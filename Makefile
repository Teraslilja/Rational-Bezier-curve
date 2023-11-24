#
# (C) Matti Lehtonen 2023
#
INCLUDE_DIRS = -I ./src/ -I ./test/

SYSTEM_HEADERS = cassert cmath cstdint functional iostream optional ostream span string_view type_traits variant vector

BASE_FILES = point_base point2 point3 bernstein_polynomials newton_raphson control_point calculation_interface validation_interface calculation validation rational_bezier
SOURCES = $(addprefix src/,$(addsuffix .cpp,$(BASE_FILES)))
OBJECTS = $(addsuffix .o,$(BASE_FILES))

HEADERS = $(addprefix src/,$(addsuffix .hpp,$(BASE_FILES)))

TESTS = test/main.cc test/point_tests.cc test/control_point_tests.cc test/bernstein_polynomials_tests.cc test/newton_raphson_tests.cc test/rational_bezier_tests.cc

GTEST_DIR = /usr/src/googletest/googletest
GTEST_SOURCE = -I $(GTEST_DIR) -I $(GTEST_DIR)/include $(GTEST_DIR)/src/gtest-all.cc
GTEST_OBJECT = gtest-all.o

LIBS=
#LIBS+=-l gtest        #  libgtest-dev
#LIBS+=-l gmock      #  libgmock-dev
#LIBS+=-l benchmark  #  libbenchmark-dev

# Which compiler?
CC=g++-13
#CC=clang++-14 -flto=full

rb_test: gtest-tool $(TESTS) $(SOURCES) $(HEADERS)
	@-rm rb_test 2> /dev/null
	$(CC) --std=c++20 -Wall -Wextra -Werror -pedantic -O3 -DNDEBUG $(INCLUDE_DIRS) -xc++ $(TESTS) $(SOURCES) $(GTEST_SOURCE) $(LIBS) -o rb_test

test: gtest-tool rb_test
	./rb_test

rb_coverage: coverage-tool $(TESTS) $(SOURCES) $(HEADERS)
	@-rm *.gcno *.gcda 2> /dev/null
	$(CC) --std=c++20 -Wall -Wextra -Werror -pedantic -Og -DNDEBUG $(INCLUDE_DIRS) -g -finline-limit=0 --coverage -fprofile-abs-path -xc++ $(TESTS) $(SOURCES) $(GTEST_SOURCE) $(LIBS) -o rb_coverage

rb_coverage.info: coverage-tool rb_coverage
	@-rm *.gcda rb_coverage.info 2> /dev/null
	./rb_coverage
	lcov --base-directory . --no-external --rc lcov_branch_coverage=1 --directory . --exclude "*.cc" --capture --output-file rb_coverage.info

rb_cleaned.info: coverage-tool rb_coverage.info
	@-rm rb.info 2> /dev/null
	lcov --rc lcov_branch_coverage=1 --add-tracefile rb_coverage.info --exclude "*.cc" --output-file rb_cleaned.info

coverage/index.html: coverage-tool rb_cleaned.info
	@-rm -rf coverage/ 2> /dev/null
	genhtml --show-details --demangle-cpp --branch-coverage --prefix $(realpath ..) rb_cleaned.info --output-directory coverage/

coverage: coverage-tool rb_coverage rb_cleaned.info coverage/index.html
	@echo ""
	@echo "See detailed coverage report at file " $(realpath ./coverage/index.html)

doc/html/index.html: doxygen-tool Doxyfile $(SOURCES) $(HEADERS)
	doxygen Doxyfile

doc: doxygen-tool doc/html/index.html $(SOURCES) $(HEADERS)
	@echo "See detailed documentation at file " $(realpath ./doc/html/index.html)

all: doc test coverage

clean:
	@-rm -rf rb_test rb_coverage *.gcno *.gcda *.info coverage/ doc/ 2> /dev/null

format: format-tool $(TESTS) $(SOURCES) $(HEADERS)
	clang-format -i $(SOURCES) $(HEADERS) $(TESTS)
	clang-format -i $(TESTS)

cppcheck: cppcheck-tool $(TESTS) $(SOURCES) $(HEADERS)
	cppcheck --platform=unix64 --language=c++ --std=c++20 $(SOURCES) $(HEADERS)
	cppcheck --platform=unix64 --language=c++ --std=c++20 $(TESTS)

.PHONY: format test coverage doc clean


# Install Google Test?
gtest-tool: /usr/local/lib/libgtest.a
	@echo "Install Google test library: sudo apt-get install libgtest-dev"

# Install Google Mock?
gmock-tool: /usr/local/lib/libgmock.a
	@echo "Install Google mock library: sudo apt-get install libgmock-dev"

# Install Google Benchmark?
benchmark-tool: /usr/lib/x86_64-linux-gnu/libbenchmark.so
	@echo "Install Google benchmark library: sudo apt-get install libbenchmark-dev"

# Install cppcheck?
cppcheck-tool: /usr/bin/cppcheck
	@echo "Install cppcheck: sudo apt-get install cppcheck"

# Install coverage tool?
coverage-tool: /usr/bin/genhtml /usr/bin/lcov
	@echo "Install lcov: sudo apt-get install lcov"

# Install doxygen tool?
doxygen-tool:
	@echo "Install doxygen: sudo apt-get install doxygen doxygen-gui"

# Install clang-format tool?
format-tool: /usr/bin/clang-format
	@echo "Install clang-format: sudo apt-get install clang-format"
