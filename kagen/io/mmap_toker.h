#pragma once

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cctype>
#include <stdexcept>
#include <string>

namespace kagen {
class MappedFileToker {
public:
    explicit MappedFileToker(const std::string& filename) {
        fd_       = OpenFile(filename);
        position_ = 0;
        length_   = FileSize(fd_);
        contents_ = static_cast<char*>(mmap(nullptr, length_, PROT_READ, MAP_PRIVATE, fd_, 0));
        if (contents_ == MAP_FAILED) {
            close(fd_);
            throw std::runtime_error("mmap failed (bad filename or empty file?)");
        }
    }

    ~MappedFileToker() {
        munmap(contents_, length_);
        close(fd_);
    }

    void Reset() {
        position_ = 0;
    }

    void Seek(const std::size_t position) {
        position_ = position;
    }

    void SkipSpaces() {
        while (ValidPosition() && (Current() == ' ' || Current() == '\t')) {
            Advance();
        }
    }

    void SkipLine() {
        while (ValidPosition() && Current() != '\n') {
            Advance();
        }
        if (ValidPosition()) {
            Advance();
        }
    }

    inline std::uint64_t ScanUnsigned() {
        std::uint64_t number = 0;
        while (TestInt()) {
            const int digit = Current() - '0';
            number          = number * 10 + digit;
            Advance();
        }
        SkipSpaces();
        return number;
    }

    void SkipInt() {
        while (TestInt()) {
            Advance();
        }
        SkipSpaces();
    }

    bool TestInt() {
        return ValidPosition() && std::isdigit(Current());
    }

    bool TestString(const char* str) {
        std::size_t pos   = position_;
        bool        match = true;
        std::size_t i     = 0;

        while (str[i] != '\0') {
            if (!ValidPosition() || str[i] != Current()) {
                match = false;
                break;
            }
            Advance();
            ++i;
        }

        position_ = pos;
        return match;
    }

    bool TestChar(const char ch) {
        return ValidPosition() && Current() == ch;
    }

    bool ConsumeChar(const char ch) {
        if (TestChar(ch)) {
            Advance();
            return true;
        }
        return false;
    }

    [[nodiscard]] bool ValidPosition() const {
        return position_ < length_;
    }

    [[nodiscard]] char Current() const {
        return contents_[position_];
    }

    void Advance() {
        ++position_;
    }

    [[nodiscard]] std::size_t Position() const {
        return position_;
    }

    [[nodiscard]] std::size_t Marked() const {
        return mark_;
    }

    void Mark() {
        mark_ = Position();
    }

    [[nodiscard]] std::size_t Length() const {
        return length_;
    }

private:
    static int OpenFile(const std::string& filename) {
        const int file = open(filename.c_str(), O_RDONLY);
        if (file < 0) {
            throw std::runtime_error("cannot open input file");
        }
        return file;
    }

    static std::size_t FileSize(const int fd) {
        struct stat file_info {};
        fstat(fd, &file_info);
        return static_cast<std::size_t>(file_info.st_size);
    }

    int         fd_;
    std::size_t position_;
    std::size_t length_;
    char*       contents_;
    std::size_t mark_;
};
} // namespace kagen
