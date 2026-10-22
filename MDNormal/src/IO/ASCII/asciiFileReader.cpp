/**
 * @file asciiFileReader.cpp
 * @brief Implementation of the asciiFileReader class.
 *
 * @author Federico Celauro
 * @date   Mar 11, 2026
 */

#include "../../include/IO/ASCII/asciiFileReader.h"

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

asciiFileReader::asciiFileReader(const std::string& path,
                                const std::string& fileName) {
    initFilePath(path, fileName);
}

asciiFileReader::~asciiFileReader() {
    if (this->in.is_open()) {
        this->in.close();
    }
}

// ---------------------------------------------------------------------------
// Stream accessors
// ---------------------------------------------------------------------------

std::ifstream& asciiFileReader::getPointer(void) {
    return in;
}

void asciiFileReader::setPointer(std::ifstream& outptr) {
    in = std::move(outptr);
}

// ---------------------------------------------------------------------------
// Path accessors
// ---------------------------------------------------------------------------

fs::path asciiFileReader::getPath(void) {
    return filePath;
}

void asciiFileReader::setPath(const fs::path& FilePath) {
    filePath = FilePath;
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

void asciiFileReader::initFilePath(const std::string& path,
                                   const std::string& fileName) {
    this->filePath = fs::path(path) / fileName;
    this->in.open(filePath);
    ensureOpen();
}

void asciiFileReader::ensureOpen() {
    if (!in.is_open()) {
        throw std::runtime_error("File not open: " + filePath.string());
    }
}

bool asciiFileReader::isFileReadable(void) const {
    if (in) {
        return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Bulk read operations
// ---------------------------------------------------------------------------

std::vector<std::string> asciiFileReader::getAllLines(void) {
    ensureOpen();

    std::vector<std::string> lineVec;
    std::string line;

    while (std::getline(in, line)) {
        lineVec.push_back(line);
    }

    resetToBeginning();
    return lineVec;
}

std::vector<std::string> asciiFileReader::getAllLines(bool skipFirstRow) {
    ensureOpen();

    std::vector<std::string> lineVec;
    std::string line;

    // Optionally skip the CSV header row
    if (skipFirstRow) {
        std::getline(in, line);
    }

    while (std::getline(in, line)) {
        lineVec.push_back(line);
    }

    resetToBeginning();
    return lineVec;
}

/**
 * @brief Reads the entire file into a single string using rdbuf().
 *
 * This is faster than line-by-line reading when the full content is needed.
 */
std::string asciiFileReader::getAllLinesOptimized(void) {
    ensureOpen();

    std::stringstream buffer;
    buffer << in.rdbuf();
    resetToBeginning();
    return buffer.str();
}

// ---------------------------------------------------------------------------
// Line-level read operations
// ---------------------------------------------------------------------------

/**
 * @brief Returns the line at a specific 0-based row index.
 *
 * If @p targetRow is before the current cursor position, the stream is
 * reset to the beginning before seeking forward.
 */
std::string asciiFileReader::getLineAt(int targetRow) {
    ensureOpen();

    if (targetRow < 0) {
        throw std::out_of_range("Negative row index is not allowed.");
    }
    if (targetRow < rowIndex) {
        resetToBeginning();
    }

    while (rowIndex < targetRow) {
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        rowIndex++;
    }

    std::string line;
    if (!std::getline(in, line)) {
        throw std::out_of_range("Desired row is beyond the end of the file.");
    }

    rowIndex = targetRow + 1;
    return line;
}

std::string asciiFileReader::getNextLine(void) {
    ensureOpen();

    std::string line;
    if (std::getline(this->in, line)) {
        rowIndex++;
        return line;
    }
    throw std::out_of_range("Desired row is beyond the end of the file.");
}

/**
 * @brief Returns multiple consecutive lines starting at @p targetRow.
 *
 * @param targetRow  0-based starting row index.
 * @param numRows    Number of rows to read.
 * @throws std::out_of_range if the requested range exceeds the file length.
 */
std::vector<std::string> asciiFileReader::getMultipleLines(int targetRow, int numRows) {
    ensureOpen();

    if (targetRow < 0) {
        throw std::out_of_range("Negative row index is not allowed.");
    }
    if (numRows <= 0) {
        throw std::out_of_range("Number of rows must be positive.");
    }

    if (targetRow < rowIndex) {
        resetToBeginning();
    }

    while (rowIndex < targetRow) {
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        rowIndex++;
    }

    std::vector<std::string> lineVec;
    lineVec.reserve(numRows);

    std::string line;
    while (lineVec.size() < static_cast<size_t>(numRows) && std::getline(in, line)) {
        lineVec.push_back(line);
    }

    if (lineVec.size() < static_cast<size_t>(numRows)) {
        throw std::out_of_range("Requested range exceeds the total number of rows in the file.");
    }

    rowIndex += lineVec.size();
    return lineVec;
}

// ---------------------------------------------------------------------------
// Cursor management
// ---------------------------------------------------------------------------

int asciiFileReader::getRowIndex() const {
    return rowIndex;
}

/**
 * @brief Resets the stream to the beginning of the file.
 *
 * Clears any EOF/error flags and moves the get-pointer to position 0.
 */
void asciiFileReader::resetToBeginning() {
    this->in.clear();                 // Clear EOF/error flags
    this->in.seekg(0, std::ios::beg); // Move cursor to the start
    this->rowIndex = 0;
}
