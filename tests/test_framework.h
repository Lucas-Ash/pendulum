#pragma once

#include <cmath>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace testfw {

class Failure : public std::runtime_error {
public:
    explicit Failure(const std::string& message)
        : std::runtime_error(message) {}
};

struct TestCase {
    std::string name;
    std::function<void()> fn;
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

inline void register_test(const std::string& name, std::function<void()> fn) {
    registry().push_back(TestCase{name, std::move(fn)});
}

template <typename T>
std::string to_string_any(const T& value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

inline std::string at(const char* file, int line) {
    std::ostringstream oss;
    oss << file << ":" << line;
    return oss.str();
}

inline std::string format_seconds(double seconds) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(seconds >= 10.0 ? 1 : 2) << seconds << "s";
    return oss.str();
}

inline int run_all() {
    const size_t total = registry().size();
    int failures = 0;
    const auto suite_start = std::chrono::steady_clock::now();

    for (size_t i = 0; i < total; ++i) {
        const auto& test = registry()[i];
        std::cout << "[RUN  " << (i + 1) << "/" << total << "] " << test.name << "\n";
        std::cout.flush();

        const auto test_start = std::chrono::steady_clock::now();
        try {
            test.fn();
            std::cout << "[PASS " << (i + 1) << "/" << total << "] " << test.name;
        } catch (const Failure& ex) {
            ++failures;
            std::cout << "[FAIL " << (i + 1) << "/" << total << "] " << test.name
                      << " - " << ex.what();
        } catch (const std::exception& ex) {
            ++failures;
            std::cout << "[FAIL " << (i + 1) << "/" << total << "] " << test.name
                      << " - unexpected exception: " << ex.what();
        } catch (...) {
            ++failures;
            std::cout << "[FAIL " << (i + 1) << "/" << total << "] " << test.name
                      << " - unknown exception";
        }

        const auto test_end = std::chrono::steady_clock::now();
        const double test_s =
            std::chrono::duration<double>(test_end - test_start).count();
        const double elapsed_s =
            std::chrono::duration<double>(test_end - suite_start).count();
        const double done = static_cast<double>(i + 1);
        const double avg_s = elapsed_s / done;
        const double eta_s = avg_s * static_cast<double>(total - (i + 1));
        const double pct = (done / static_cast<double>(total)) * 100.0;

        std::cout << " | test=" << format_seconds(test_s)
                  << " | progress=" << std::fixed << std::setprecision(1)
                  << pct << "%"
                  << " | elapsed=" << format_seconds(elapsed_s)
                  << " | eta~" << format_seconds(eta_s) << "\n";
    }

    std::cout << "\n";
    std::cout << "Total tests: " << total << "\n";
    std::cout << "Failures: " << failures << "\n";
    return failures == 0 ? 0 : 1;
}

}  // namespace testfw

#define TEST(NAME)                                                              \
    void NAME();                                                                \
    namespace {                                                                 \
    struct NAME##_registrar {                                                   \
        NAME##_registrar() { testfw::register_test(#NAME, NAME); }             \
    };                                                                          \
    static NAME##_registrar NAME##_registrar_instance;                          \
    }                                                                           \
    void NAME()

#define EXPECT_TRUE(EXPR)                                                       \
    do {                                                                        \
        if (!(EXPR)) {                                                          \
            throw testfw::Failure(testfw::at(__FILE__, __LINE__) +             \
                                  " EXPECT_TRUE failed: " #EXPR);              \
        }                                                                       \
    } while (0)

#define EXPECT_FALSE(EXPR) EXPECT_TRUE(!(EXPR))

#define EXPECT_EQ(A, B)                                                         \
    do {                                                                        \
        const auto _a = (A);                                                    \
        const auto _b = (B);                                                    \
        if (!(_a == _b)) {                                                      \
            throw testfw::Failure(testfw::at(__FILE__, __LINE__) +             \
                                  " EXPECT_EQ failed: got [" +                 \
                                  testfw::to_string_any(_a) + "] expected [" + \
                                  testfw::to_string_any(_b) + "]");            \
        }                                                                       \
    } while (0)

#define EXPECT_NEAR(A, B, TOL)                                                  \
    do {                                                                        \
        const double _a = static_cast<double>(A);                               \
        const double _b = static_cast<double>(B);                               \
        const double _tol = static_cast<double>(TOL);                           \
        if (std::fabs(_a - _b) > _tol) {                                        \
            throw testfw::Failure(testfw::at(__FILE__, __LINE__) +             \
                                  " EXPECT_NEAR failed: |" +                   \
                                  testfw::to_string_any(_a) + " - " +         \
                                  testfw::to_string_any(_b) + "| > " +        \
                                  testfw::to_string_any(_tol));                \
        }                                                                       \
    } while (0)

#define EXPECT_FINITE(VALUE)                                                    \
    do {                                                                        \
        const double _v = static_cast<double>(VALUE);                           \
        if (!std::isfinite(_v)) {                                               \
            throw testfw::Failure(testfw::at(__FILE__, __LINE__) +             \
                                  " EXPECT_FINITE failed for " #VALUE);        \
        }                                                                       \
    } while (0)

#define EXPECT_THROW(STMT)                                                      \
    do {                                                                        \
        bool _thrown = false;                                                   \
        try {                                                                   \
            STMT;                                                               \
        } catch (...) {                                                         \
            _thrown = true;                                                     \
        }                                                                       \
        if (!_thrown) {                                                         \
            throw testfw::Failure(testfw::at(__FILE__, __LINE__) +             \
                                  " EXPECT_THROW failed: " #STMT);             \
        }                                                                       \
    } while (0)
