//
// (C) Matti Lehtonen 2023
//

#include <gtest/gtest.h>

int
main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  int const status = RUN_ALL_TESTS();
  if (status != 0) {
    std::exit(status);
  }

  std::exit(EXIT_SUCCESS);
}
