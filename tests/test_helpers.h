#pragma once

#include <filesystem>
#include <string>
#include <vector>

struct CsvData {
    std::string header;
    std::vector<std::vector<double>> rows;
};

struct DatData {
    std::string header;
    std::vector<std::vector<double>> rows;
};

class TempDir {
public:
    TempDir();
    ~TempDir();

    const std::filesystem::path& path() const;
    std::filesystem::path child(const std::string& relative) const;
    std::filesystem::path write_file(const std::string& relative, const std::string& content) const;

private:
    std::filesystem::path path_;
};

class EnvVarGuard {
public:
    explicit EnvVarGuard(const std::string& key);
    ~EnvVarGuard();

    void set(const std::string& value) const;
    void unset() const;

private:
    std::string key_;
    bool had_original_ = false;
    std::string original_;
};

class WorkingDirGuard {
public:
    explicit WorkingDirGuard(const std::filesystem::path& target);
    ~WorkingDirGuard();

private:
    std::filesystem::path original_;
};

CsvData read_csv(const std::filesystem::path& path);
DatData read_dat(const std::filesystem::path& path);
double max_abs(const std::vector<double>& values);
double mean_abs(const std::vector<double>& values);
