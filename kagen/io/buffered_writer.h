#pragma once

#include "kagen/io.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdlib>
#include <string>

namespace kagen {
struct CreateTag {};
struct AppendTag {};

namespace tag {
constexpr CreateTag create;
constexpr AppendTag append;
} // namespace tag

template <std::size_t kBufferSize = 1024 * 1024, std::size_t kBufferSizeLimit = kBufferSize - 1024>
class BufferedTextOutput {
public:
    BufferedTextOutput(CreateTag, const std::string& filename)
        : fd_{open(filename.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR)} {
        if (fd_ < 0) {
            throw IOError("cannot write to " + filename);
        }
    }

    BufferedTextOutput(AppendTag, const std::string& filename) : fd_{open(filename.c_str(), O_WRONLY | O_APPEND)} {
        if (fd_ < 0) {
            throw IOError("cannot write to " + filename + " (this is most likely a bug)");
        }
    }

    ~BufferedTextOutput() {
        ForceFlush();
        close(fd_);
    }

    BufferedTextOutput& WriteString(const char* str) {
        for (const char* ch = str; *ch; ++ch) {
            WriteChar(*ch);
        }
        return *this;
    }

    BufferedTextOutput& WriteChar(const char ch) {
        *(buffer_pos_)++ = ch;
        return *this;
    }

    template <typename Int>
    BufferedTextOutput& WriteInt(Int value) {
        static char rev_buffer[80];

        int pos = 0;
        do {
            rev_buffer[pos++] = value % 10;
            value /= 10;
        } while (value > 0);

        while (pos > 0) {
            *(buffer_pos_++) = '0' + rev_buffer[--pos];
        }
        return *this;
    }

    BufferedTextOutput& WriteFloat(const double value) {
        int written = std::sprintf(buffer_pos_, "%.5lf", value);
        buffer_pos_ += written;
        return *this;
    }

    BufferedTextOutput& Flush() {
        if (static_cast<std::size_t>(buffer_pos_ - buffer_) >= kBufferSizeLimit) {
            ForceFlush();
        }
        return *this;
    }

private:
    void ForceFlush() {
        [[maybe_unused]] auto nbytes_writte = write(fd_, buffer_, buffer_pos_ - buffer_);
        buffer_pos_                         = buffer_;
    }

    int   fd_;
    char  buffer_[kBufferSize]{0};
    char* buffer_pos_{buffer_};
};
} // namespace kagen
