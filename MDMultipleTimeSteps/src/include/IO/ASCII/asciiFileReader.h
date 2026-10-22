/**
 * @file asciiFileReader.h
 * @brief Provides file-reading utilities for ASCII/CSV data files.
 *
 * Wraps std::ifstream with convenience methods for reading individual
 * lines, ranges of lines, or entire files.  Supports random-access
 * seeking and automatic stream management.
 *
 * @author Federico Celauro
 * @date   Mar 11, 2026
 */

#ifndef ASCIIFILEREADER_H_
#define ASCIIFILEREADER_H_

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <experimental/filesystem>
#include <limits>
#include <utility>

namespace fs = std::experimental::filesystem;

/**
 * @class asciiFileReader
 * @brief RAII wrapper around std::ifstream for sequential and random-access
 *        reading of text files.
 */
class asciiFileReader {
public:
    /// Constructs the reader and opens the file at @p path / @p fileName.
    asciiFileReader(const std::string& path, const std::string& fileName);

    virtual ~asciiFileReader();

    /// @name Stream accessors
    /// @{
    std::ifstream& getPointer(void);
    void setPointer(std::ifstream& stream);
    /// @}

    /// @name Path accessors
    /// @{
    fs::path getPath(void);
    void setPath(const fs::path& filePath);
    /// @}

    /// @name Bulk read operations
    /// @{

    /// Reads and returns all lines in the file.
    std::vector<std::string> getAllLines(void);

    /// Reads all lines, optionally skipping the first row (header).
    std::vector<std::string> getAllLines(bool skipFirstRow);

    /// Reads the entire file content into a single string (fast bulk read).
    std::string getAllLinesOptimized(void);
    /// @}

    /// @name Line-level read operations
    /// @{

    /// Returns the line at a specific 0-based row index.
    std::string getLineAt(int targetRow);

    /// Returns the next unread line and advances the internal cursor.
    std::string getNextLine(void);

    /// Returns @p numRows consecutive lines starting from @p targetRow.
    std::vector<std::string> getMultipleLines(int targetRow, int numRows);
    /// @}

    /// Returns the 0-based index of the next unread row.
    int getRowIndex() const;

    /// Returns true if the underlying stream is in a good state.
    bool isFileReadable() const;

    /// Resets the stream cursor to the beginning of the file.
    void resetToBeginning(void);

private:
    fs::path filePath;       ///< Full path to the file
    std::ifstream in;        ///< Input file stream
    int rowIndex = 0;        ///< Current row cursor position

    /// Builds the full path and opens the file.
    void initFilePath(const std::string& path, const std::string& fileName);

    /// Throws std::runtime_error if the file is not open.
    void ensureOpen(void);
};

#endif /* ASCIIFILEREADER_H_ */
