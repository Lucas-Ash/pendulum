#include "tests/test_framework.h"

#include <climits>

#include "modules/config_utils.h"
#include "tests/test_helpers.h"

TEST(ConfigUtilsTrimHandlesEdgeCases) {
    EXPECT_EQ(config_utils::trim(""), "");
    EXPECT_EQ(config_utils::trim("   \t \n"), "");
    EXPECT_EQ(config_utils::trim("  hello"), "hello");
    EXPECT_EQ(config_utils::trim("hello  "), "hello");
    EXPECT_EQ(config_utils::trim("  hello  "), "hello");
    EXPECT_EQ(config_utils::trim("a b c"), "a b c");
}

TEST(ConfigUtilsToLowerHandlesEdgeCases) {
    EXPECT_EQ(config_utils::to_lower("MiXeD"), "mixed");
    EXPECT_EQ(config_utils::to_lower("lower"), "lower");
    EXPECT_EQ(config_utils::to_lower(""), "");
    EXPECT_EQ(config_utils::to_lower("UPPER"), "upper");
}

TEST(ConfigUtilsParseDoubleValidAndInvalid) {
    double value = 0.0;
    EXPECT_TRUE(config_utils::parse_double("42", value));
    EXPECT_NEAR(value, 42.0, 1e-12);

    EXPECT_TRUE(config_utils::parse_double("-3.125", value));
    EXPECT_NEAR(value, -3.125, 1e-12);

    EXPECT_TRUE(config_utils::parse_double("1e-6", value));
    EXPECT_NEAR(value, 1e-6, 1e-18);

    EXPECT_FALSE(config_utils::parse_double("abc", value));
    EXPECT_FALSE(config_utils::parse_double("12abc", value));
    EXPECT_FALSE(config_utils::parse_double("", value));
    EXPECT_FALSE(config_utils::parse_double("1e309", value));
}

TEST(ConfigUtilsParseIntValidAndInvalid) {
    int value = 0;
    EXPECT_TRUE(config_utils::parse_int("123", value));
    EXPECT_EQ(value, 123);

    EXPECT_TRUE(config_utils::parse_int("-456", value));
    EXPECT_EQ(value, -456);

    EXPECT_TRUE(config_utils::parse_int(std::to_string(INT_MAX), value));
    EXPECT_EQ(value, INT_MAX);

    EXPECT_FALSE(config_utils::parse_int("3.14", value));
    EXPECT_FALSE(config_utils::parse_int("abc", value));
    EXPECT_FALSE(config_utils::parse_int("", value));
}

TEST(ConfigUtilsParseBoolSupportedForms) {
    bool out = false;

    EXPECT_TRUE(config_utils::parse_bool("true", out));
    EXPECT_TRUE(out);
    EXPECT_TRUE(config_utils::parse_bool("YES", out));
    EXPECT_TRUE(out);
    EXPECT_TRUE(config_utils::parse_bool("On", out));
    EXPECT_TRUE(out);
    EXPECT_TRUE(config_utils::parse_bool("1", out));
    EXPECT_TRUE(out);

    EXPECT_TRUE(config_utils::parse_bool("false", out));
    EXPECT_FALSE(out);
    EXPECT_TRUE(config_utils::parse_bool("No", out));
    EXPECT_FALSE(out);
    EXPECT_TRUE(config_utils::parse_bool("off", out));
    EXPECT_FALSE(out);
    EXPECT_TRUE(config_utils::parse_bool("0", out));
    EXPECT_FALSE(out);

    EXPECT_FALSE(config_utils::parse_bool("maybe", out));
}

TEST(ConfigUtilsParseStringHandlesQuotedAndUnquoted) {
    EXPECT_EQ(config_utils::parse_string("value"), "value");
    EXPECT_EQ(config_utils::parse_string("\"quoted value\""), "quoted value");
    EXPECT_EQ(config_utils::parse_string("'quoted value'"), "quoted value");
    EXPECT_EQ(config_utils::parse_string("   \"trim me\"   "), "trim me");
    EXPECT_EQ(config_utils::parse_string(""), "");
}

TEST(ConfigUtilsResolveOutputPathWithQaTestEnv) {
    EnvVarGuard guard("QA_TEST");
    guard.set("1");

    const std::string original = "some/dir/file.csv";
    EXPECT_EQ(config_utils::resolve_output_path(original), original);
}

TEST(ConfigUtilsResolveOutputPathWithoutQaTestEnv) {
    EnvVarGuard guard("QA_TEST");
    guard.unset();

    TempDir temp;
    WorkingDirGuard wd(temp.path());

    const std::string resolved = config_utils::resolve_output_path("x/y/data.csv");
    EXPECT_EQ(resolved, "outputs/data.csv");
    EXPECT_TRUE(std::filesystem::exists(temp.child("outputs")));
}
