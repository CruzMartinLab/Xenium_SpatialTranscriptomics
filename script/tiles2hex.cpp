#include "tiles2bins.hpp"

int32_t cmdTiles2HexTxt(int32_t argc, char** argv) {

    std::string inTsv, inIndex, outFile, tmpDir, dictFile;
    std::vector<std::string> anchorFiles;
    std::vector<float> radius;
    int nThreads = 1, debug = 0, verbose = 1000000;
    int icol_x, icol_y, icol_feature;
    double hexSize;
    std::vector<int32_t> icol_ints;
    std::vector<int32_t> min_counts;
    bool noBackground = false;

    ParamList pl;
    // Input Options
    pl.add_option("in-tsv", "Input TSV file. Header must begin with #", inTsv)
        .add_option("in-index", "Input index file", inIndex)
        .add_option("icol-x", "Column index for x coordinate (0-based)", icol_x)
        .add_option("icol-y", "Column index for y coordinate (0-based)", icol_y)
        .add_option("icol-feature", "Column index for feature (0-based)", icol_feature)
        .add_option("feature-dict", "If feature column is not integer, provide a dictionary/list of all possible values", dictFile)
        .add_option("icol-int", "Column index for integer values (0-based)", icol_ints)
        .add_option("anchor-files", "Anchor files", anchorFiles)
        .add_option("radius", "Radius for each set of anchors", radius)
        .add_option("hex-size", "Hexagon size (size length)", hexSize)
        .add_option("temp-dir", "Directory to store temporary files", tmpDir)
        .add_option("threads", "Number of threads to use (default: 1)", nThreads);
    // Output Options
    pl.add_option("out", "Output TSV file", outFile)
        .add_option("min-count", "Minimum count for each integer column, applied with OR", min_counts)
        .add_option("ignore-background", "Ignore pixels not within radius of any of the anchors", noBackground)
        .add_option("verbose", "Verbose", verbose)
        .add_option("debug", "Debug", debug);

    try {
        pl.readArgs(argc, argv);
        pl.print_options();
    } catch (const std::exception &ex) {
        std::cerr << "Error parsing options: " << ex.what() << "\n";
        pl.print_help();
        return 1;
    }

    HexGrid hexGrid(hexSize);
    TileReader tileReader(inTsv, inIndex);
    if (!tileReader.isValid()) {
        error("Error opening input file: %s", inTsv.c_str());
        return 1;
    }
    lineParser parser(icol_x, icol_y, icol_feature, icol_ints, dictFile);
    if (parser.n_ints == 0) {
        error("No integer columns specified");
    }

    if (anchorFiles.empty()) {
        Tiles2Hex tiles2Hex(nThreads, tmpDir, outFile, hexGrid, tileReader, parser, min_counts);
        if (!tiles2Hex.run()) {
            return 1;
        }
        tiles2Hex.writeMetadata();
    } else {
        Tiles2UnitsByAnchor tiles2Hex(nThreads, tmpDir, outFile, hexGrid, tileReader, parser, anchorFiles, radius, min_counts, noBackground);
        if (!tiles2Hex.run()) {
            return 1;
        }
        tiles2Hex.writeMetadata();
    }

    notice("Processing completed. Output is written to %s", outFile.c_str());

    return 0;
}
