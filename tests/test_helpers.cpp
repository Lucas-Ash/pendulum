#include "tests/test_helpers.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <unistd.h>

namespace {

std::string unique_suffix() {
    const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::ostringstream oss;
    oss << ::getpid() << "_" << now;
    return oss.str();
}

}  // namespace

TempDir::TempDir() {
    path_ = std::filesystem::temp_directory_path() / ("pendulum_tests_" + unique_suffix());
    std::filesystem::create_directories(path_);
}

TempDir::~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(path_, ec);
}

const std::filesystem::path& TempDir::path() const {
    return path_;
}

std::filesystem::path TempDir::child(const std::string& relative) const {
    return path_ / relative;
}

std::filesystem::path TempDir::write_file(const std::string& relative, const std::string& content) const {
    const std::filesystem::path out = child(relative);
    std::filesystem::create_directories(out.parent_path());
    std::ofstream file(out);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to write file: " + out.string());
    }
    file << content;
    return out;
}

EnvVarGuard::EnvVarGuard(const std::string& key)
    : key_(key) {
    const char* existing = std::getenv(key_.c_str());
    if (existing != nullptr) {
        had_original_ = true;
        original_ = existing;
    }
}

EnvVarGuard::~EnvVarGuard() {
    if (had_original_) {
        ::setenv(key_.c_str(), original_.c_str(), 1);
    } else {
        ::unsetenv(key_.c_str());
    }
}

void EnvVarGuard::set(const std::string& value) const {
    ::setenv(key_.c_str(), value.c_str(), 1);
}

void EnvVarGuard::unset() const {
    ::unsetenv(key_.c_str());
}

WorkingDirGuard::WorkingDirGuard(const std::filesystem::path& target)
    : original_(std::filesystem::current_path()) {
    std::filesystem::current_path(target);
}

WorkingDirGuard::~WorkingDirGuard() {
    std::error_code ec;
    std::filesystem::current_path(original_, ec);
}

CsvData read_csv(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open CSV: " + path.string());
    }

    CsvData data;
    std::getline(file, data.header);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<double> row;
        std::stringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) {
            row.push_back(std::stod(token));
        }
        data.rows.push_back(std::move(row));
    }

    return data;
}

DatData read_dat(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open DAT: " + path.string());
    }

    DatData data;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            data.header = line;
            continue;
        }

        std::vector<double> row;
        std::stringstream ss(line);
        double value = 0.0;
        while (ss >> value) {
            row.push_back(value);
        }
        if (!row.empty()) {
            data.rows.push_back(std::move(row));
        }
    }
    return data;
}

double max_abs(const std::vector<double>& values) {
    double out = 0.0;
    for (double v : values) {
        out = std::max(out, std::fabs(v));
    }
    return out;
}

double mean_abs(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    double acc = 0.0;
    for (double v : values) {
        acc += std::fabs(v);
    }
    return acc / static_cast<double>(values.size());
}
