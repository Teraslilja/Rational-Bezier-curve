//
// (C) Matti Lehtonen 2023
//

/**
 *  @file validation_interface.cpp This file contains implementation related to delegation interface
 */

#include "validation_interface.hpp"

namespace curve::bezier::rational {

// LCOV_EXCL_START
static inline constexpr std::optional<std::string_view const> toString(ValidityIssue const data) {
  constexpr std::string_view const issue1{"ISSUE_U_IS_INVALID"};
  constexpr std::string_view const issue2{"ISSUE_NOT_ENOUGHT_CONTROL_POINTS"};
  constexpr std::string_view const issue3{"ISSUE_BAD_COMBINATION_OF_WEIGHTS"};
  constexpr std::string_view const issue4{"ISSUE_OUT_OF_HEAP_MEMORY"};
  constexpr std::string_view const issue5{"ISSUE_BAD_CONTROLPOINT_WEIGHT"};
  constexpr std::string_view const issue6{"ISSUE_BAD_POINT"};

  switch (data) {
  case ValidityIssue::ISSUE_U_IS_INVALID:
    return std::make_optional(issue1);
  case ValidityIssue::ISSUE_NOT_ENOUGHT_CONTROL_POINTS:
    return std::make_optional(issue2);
  case ValidityIssue::ISSUE_BAD_COMBINATION_OF_WEIGHTS:
    return std::make_optional(issue3);
  case ValidityIssue::ISSUE_OUT_OF_HEAP_MEMORY:
    return std::make_optional(issue4);
  case ValidityIssue::ISSUE_BAD_CONTROLPOINT_WEIGHT:
    return std::make_optional(issue5);
  case ValidityIssue::ISSUE_BAD_POINT:
    return std::make_optional(issue6);
  }
  return std::nullopt;
}

std::ostream &operator<<(std::ostream &out, ValidityIssue const data) {
  auto const str = toString(data);
  if (str.has_value()) {
    out << str.value().substr();
  } else {
    out << "'invalid issue <" << static_cast<std::uint16_t>(data) << ">'";
  }
  return out;
}

std::ostream &operator<<(std::ostream &out, std::optional<ValidityIssue> const &data) {
  if (data.has_value()) {
    out << data.value();
  } else {
    out << "'no value'";
  }
  return out;
}
// LCOV_EXCL_STOP

} // namespace curve::bezier::rational
