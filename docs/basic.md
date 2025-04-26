# Pixel Level Factor Analysis Example

This example demonstrates how to perform pixel level factor analysis with punkst, achieving similar results to FICTURE (2024) with improved efficiency.

We will explain how to generate a full workflow using a template Makefile, then explains the steps taken inside the workflow.

## Generic input format and example data

There is a small example data in `punkst/examples/data`. Decompress the the data to your local directory `tar xvzf example_data.tar.gz -C /path/to/parent/dir`.

See [Input](./input/index.md) for details on starting from raw data from different platforms.

We need three pieces of information:

- The main input file is a TSV file with X coordinate, Y coordinate, feature, and count columns. (`example_data/transcripts.tsv`) It can have other columns (which will be ignored), and it may or may not have headers (see Step 1). We will only use this file directly in step 1.

- A file containing a list of genes to use, one per line. If your file contains multiple columns, the first column delimited by any whitespace will be used. (`example_data/features.txt`)

- The coordinate range: xmin, xmax, ymin, ymax. (For the example data you can find them in `example_data/coordinate_minmax.tsv`)



## Use the example Makefile template

We provide a template Makefile and config file in `punkst/examples` to generate the full workflow of FICTURE.

You can copy `punkst/examples/basic/generic/config.json` to your own directory modify the data path and parameters, then use `punkst/ext/py/generate_workflow.py` to generate a data-specific Makefile for your task.

The python script also generates a bash script that can be submitted as a slurm job. If you are not using slurm just ignore the parameters in the "job"  section of the config and run the generation script without the `-o` option.

```bash
# set repopath to the path of the punkst repo
python ${repopath}/ext/py/generate_workflow.py \
  -c config.json -o run.sh -m Makefile \
  -t ${repopath}/examples/basic/generic/Makefile
```

You can check the generated workflow before execution by
```bash
make -f Makefile --dry-run
```

Then `make -f Makefile` exectutes the workflow.

### Parameters in config.json

(The parameters in the example config file works for the example data, where the coordinates are in microns.)

- `"datadir"`: the path to store all output

- `"tmpdir"`: the path to store temporary files (those files will be deleted automatically by the program). This directory must be empty or creatable.

- `"transcripts"`: a tsv file with X coordinate, Y coordinate, gene/transcript name, and count columns. There could be other columns but they will be ignored.

- Specify the 0-based column indices in "transcripts" for X coordinate, Y coordinate, feature, and count: `"icol_x"`, `"icol_y"`, `"icol_feature"`, and `"icol_count"`. If the input file contains headers, set "skip" to the number of lines to skip.

- `"features"`: a file that contains a list of genes to use. Rows in the transcripts file where the gene name on the column specified by `icol_feature` does not match any names in this file will be ignored. (If your feature list contains multiple columns, the first column delimited by any whitespace will be used.)

- `"tilesize"`: we store and process data by square tiles, this parameter specifies the size length of the tiles in the same unit as your coordinates. Tile sizes affect the memory usage and (perhaps less so) run time, we've been using 500$\mu$m for all of our experiments.

- `"hexgrids"` (list): this is center-to-center distance of the hexagonal grid used for training the model. The best value depends on your data density. We've been using 12~18$\mu$m, but you might want to use a larger value if your data has very low transcript density.


- `"topics"` (list): the number of topics (factors) to learn.

- `"pixhex"`: often set to be the same as "hexgrids" or slightly smaller.

- `"nmove"`: "pixhex" divided by "nmove" is the distance between adjacent anchor points in the algorithm. We recommend pixhex/nmove to be around $4~6\mu m$ for high resolution results.

- `"res"`: the resolution for pixel level inference (pixels within this distance will be grouped together in inference). We've been using $0.5\mu m$.

- `"scale"`: this only controls the visualization of pixel level results. The coordinate values divided by scale will be the "pixel" indices in the image. If your coordinates are in microns and you want $0.5 \mu m$ to be one pixel in the image, set scale to 0.5. For Visium HD where the data resolution is $2 \mu m$, you probably want to set scale to 2.

- `"xmin", "xmax", "ymin", "ymax"`: the range of coordinates.

- Section `"job"`: only for slurm users. Those are just slurm job parameters to create a job script to wrap aroun the Makefile. You probably don't need this, just for convenience. You can include additional commands by setting "extra_lines".


## Step by step

### Setup

First, set up the environment variables:

```bash
threads=4 # Number of threads for parallel processing
tmpdir=/path/to/tmp # Directory for temporary files (must be empty or creatable)
path=/path/to/your_data # Path to your data directory
```

### Step 1: Group pixels to tiles

Group pixels into non-overlapping square tiles for faster processing:

```bash
punkst pts2tiles --in-tsv ${path}/transcripts.tsv \
  --icol-x 0 --icol-y 1 --skip 0 \
  --tile-size 500 \
  --temp-dir ${tmpdir} --threads ${threads} \
  --out-prefix ${path}/transcripts.tiled
```

Key parameters:
- `--icol-x`, `--icol-y`: Column indices for X and Y coordinates (0-based)
- If your input file has a header, use `--skip 1` to skip the first (or more) lines
- `--tile-size`: Size (side length) of the square tiles

[Detailed documentation for pts2tiles](./modules/pts2tiles.md)

### Step 2: Create hexagonal units

Group pixels into non-overlapping hexagons:

```bash
punkst tiles2hex --in-tsv ${path}/transcripts.tiled.tsv \
  --in-index ${path}/transcripts.tiled.index \
  --feature-dict ${path}/features.txt \
  --icol-x 0 --icol-y 1 --icol-feature 2 --icol-int 3 \
  --min-count 20 --hex-size 7 \
  --out ${path}/hex.txt \
  --temp-dir ${tmpdir} --threads ${threads}
```

Key parameters:
- `--icol-feature`, `--icol-int`: Column indices for feature and count(s)
- `--hex-size`: Side length of the hexagons
- `--min-count`: Minimum count for a hexagon to be included

Shuffle the output for training:
```bash
sort -k1,1 --parallel ${threads} ${path}/hex.txt > ${path}/hex.randomized.txt
rm ${path}/hex.txt
```

[Detailed documentation for tiles2hex](./modules/tiles2hex.md)

### Step 3: Run LDA on hexagon data

Perform Latent Dirichlet Allocation on the hexagon data:

```bash
punkst lda4hex --in-data ${path}/hex.randomized.txt \
  --in-meta ${path}/hex.json \
  --n-topics 12 \
  --n-epochs 2 --min-count-train 50 \
  --out-prefix ${path}/hex.lda --transform \
  --threads ${threads} --seed 1
```

Key parameters:
- `--n-topics`: Number of topics (factors) to learn
- `--transform`: Generate transform results after model fitting

[Detailed documentation for lda4hex](./modules/lda4hex.md)

### Step 4: Decode pixels with the model

Annotate each pixel with top factors and their probabilities:

```bash
punkst pixel-decode --model ${path}/hex.lda.model.tsv \
  --in-tsv ${path}/transcripts.tiled.tsv \
  --in-index ${path}/transcripts.tiled.index \
  --icol-x 0 --icol-y 1 --icol-feature 2 --icol-val 3 \
  --hex-grid-dist 12 --n-moves 2 \
  --pixel-res 0.5 \
  --out-pref ${path}/pixel.decode \
  --temp-dir ${tmpdir} \
  --threads ${threads} --seed 1 --output-original
```

Key parameters:
- `--model`: Model file created by lda4hex
- `--hex-grid-dist`: Center-to-center distance of the hexagonal grid
- `--n-moves`: Number of sliding moves to generate anchors
- `--pixel-res`: Resolution for the analysis (in the same unit as coordinates)
- `--output-original`: Write each transcript/input pixel as a separate line in the output. This will be slower and generates a bigger file, so only use it if matching the inference with the original input is useful. (Excluding this flag for Visium HD data is more sensible)

[Detailed documentation for pixel-decode](./modules/pixel-decode.md)

### Step 5: Visualize the results

Visualize the pixel decoding results:

Optional: choose a color table based on the intermediate results. Otherwise, you need to create a RGB table with the following columns: R, G, B (including the header). So each row represents the RGB color for a factor, with integer values from 0 to 255. (Python dependency: [jinja2](https://pypi.org/project/Jinja2/), pandas, matplotlib.)

```bash
python punkst/ext/py/color_helper.py --input ${path}/hex.lda.results.tsv --output ${path}/color
```

Generate an image for the pixel level factor assignment
```bash
punkst draw-pixel-factors --in-tsv ${path}/pixel.decode.tsv \
  --in-color ${path}/color.rgb.tsv \
  --out ${path}/pixel.png \
  --scale 1 \
  --xmin ${xmin} --xmax ${xmax} --ymin ${ymin} --ymax ${ymax}
```

Key parameters:
- `--in-color`: TSV file with RGB colors for each factor
- `--scale`: Scales input coordinates to pixels in the output image (2 means 2 coordinate units = 1 pixel)
- `--xmin`, `--xmax`, `--ymin`, `--ymax`: Range of coordinates to visualize

[Detailed documentation for visualization](./modules/visualization.md)
