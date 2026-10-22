/**
 * @file asciiFileWriter.h
 * @brief Provides file-writing utilities for ASCII/CSV data files.
 *
 * Wraps std::ofstream with convenience methods for writing single lines,
 * vectors of lines, and switching between overwrite and append modes.
 *
 * @author Federico Celauro
 * @date   Mar 12, 2026
 */

#ifndef ASCIIFILEWRITER_H_
#define ASCIIFILEWRITER_H_

#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

/**
 * @class asciiFileWriter
 * @brief RAII wrapper around std::ofstream for writing text files
 *        in either overwrite or append mode.
 */
class asciiFileWriter {
public:
    /// @name Constructors & Destructor
    /// @{

    /// Opens the file in overwrite mode.
    asciiFileWriter(const std::string& path, const std::string& fileName);

    /// Opens the file with a specific open mode.
    asciiFileWriter(const std::string& path, const std::string& fileName,
                    const std::ios_base::openmode mode);

    /// Opens the file in append mode if @p append is true, overwrite otherwise.
    asciiFileWriter(const std::string& path, const std::string& fileName,
                    const bool append);

    virtual ~asciiFileWriter();
    /// @}

    /// @name Stream accessors
    /// @{
    std::ofstream& getPointer(void);
    void setPointer(std::ofstream& stream);
    /// @}

    /// @name Path accessors
    /// @{
    fs::path getPath(void);
    void setPath(const fs::path& filePath);
    /// @}

    /// @name Write operations
    /// @{

    /// Writes a single line followed by a newline.
    void writeOnFile(const std::string& data);

    /// Writes multiple lines, each followed by a newline.
    void writeOnFile(const std::vector<std::string>& lines);
    /// @}

    /// @name Mode switching
    /// @{
    void setAppendMode(void);
    void setOverwriteMode(void);
    /// @}

private:
    std::ofstream out;   ///< Output file stream
    fs::path filePath;   ///< Full path to the output file

    /// Opens the file in overwrite mode and stores the path.
    void initFilePathOverwrite(const std::string& path, const std::string& fileName);

    /// Opens the file in append mode and stores the path.
    void initFilePathAppend(const std::string& path, const std::string& fileName);

    /// Closes any open stream and re-opens with the given mode.
    void setMode(std::ios_base::openmode mode);
};

#endif /* ASCIIFILEWRITER_H_ */
