/**
 * @file asciiFileWriter.cpp
 * @brief Implementation of the asciiFileWriter class.
 *
 * @author Federico Celauro
 * @date   Mar 12, 2026
 */

#include "../../include/IO/ASCII/asciiFileWriter.h"

// ---------------------------------------------------------------------------
// Constructors & Destructor
// ---------------------------------------------------------------------------

asciiFileWriter::asciiFileWriter(const std::string& path,
                                const std::string& fileName) {
    initFilePathOverwrite(path, fileName);
}

asciiFileWriter::asciiFileWriter(const std::string& path,
                                const std::string& fileName,
                                const std::ios_base::openmode mode) {
    filePath = fs::path(path) / fileName;
    setMode(mode);
}

asciiFileWriter::asciiFileWriter(const std::string& path,
                                const std::string& fileName,
                                const bool append) {
    if (append) {
        initFilePathAppend(path, fileName);
    } else {
        initFilePathOverwrite(path, fileName);
    }
}

asciiFileWriter::~asciiFileWriter() {
    if (this->out.is_open()) {
        try { this->out.close(); } catch (const std::ios_base::failure&) {}
    }
}

// ---------------------------------------------------------------------------
// Stream accessors
// ---------------------------------------------------------------------------

std::ofstream& asciiFileWriter::getPointer(void) {
    return out;
}

void asciiFileWriter::setPointer(std::ofstream& outptr) {
    out = std::move(outptr);
}

// ---------------------------------------------------------------------------
// Path accessors
// ---------------------------------------------------------------------------

fs::path asciiFileWriter::getPath(void) {
    return filePath;
}

void asciiFileWriter::setPath(const fs::path& FilePath) {
    filePath = FilePath;
}

// ---------------------------------------------------------------------------
// Write operations
// ---------------------------------------------------------------------------

void asciiFileWriter::writeOnFile(const std::string& data) {
    if (!(out << data << "\n")) {
        throw std::runtime_error("Failed to write to file: " + this->filePath.string());
    }
}

void asciiFileWriter::writeOnFile(const std::vector<std::string>& lineVec) {
    if (!this->out.is_open()) {
        throw std::runtime_error("Failed to write to file: " + this->filePath.string());
    }
    for (const std::string& line : lineVec) {
        this->out << line << "\n";
    }
    this->out.flush();
}

// ---------------------------------------------------------------------------
// Initialization & mode switching
// ---------------------------------------------------------------------------

void asciiFileWriter::initFilePathOverwrite(const std::string& path,
                                            const std::string& fileName) {
    this->filePath = fs::path(path) / fileName;
    setOverwriteMode();
}

void asciiFileWriter::initFilePathAppend(const std::string& path,
                                         const std::string& fileName) {
    this->filePath = fs::path(path) / fileName;
    setAppendMode();
}

void asciiFileWriter::setAppendMode(void) {
    setMode(std::ios::app);
}

void asciiFileWriter::setOverwriteMode(void) {
    setMode(std::ios::out);
}

/**
 * @brief Closes any open stream and re-opens with the specified mode.
 * @throws std::runtime_error if the file cannot be opened.
 */
void asciiFileWriter::setMode(std::ios_base::openmode mode) {
    if (this->out.is_open()) {
        this->out.close();
    }
    this->out.clear();
    this->out.exceptions(std::ios::failbit | std::ios::badbit);
    this->out.open(this->filePath, mode);

    if (!out) {
        throw std::runtime_error("Could not open file: " + this->filePath.string());
    }
}
