#
# (C) Matti Lehtonen 2023
#
INCLUDE_DIRS = -I ./src/ -I ./test/
HEADERS = src/bernstein_polynomials.hpp src/control_point.hpp src/point.hpp src/newton_raphson.hpp src/rational_bezier.hpp
SOURCES = src/point.cpp src/control_point.cpp src/bernstein_polynomials.cpp src/newton_raphson.cpp src/rational_bezier.cpp
TESTS = test/point_tests.cc test/control_point_tests.cc test/bernstein_polynomials_tests.cc test/rational_bezier_tests.cc test/newton_raphson_tests.cc test/main.cc

LIBS=-l gtest        #  libgtest-dev
#LIBS+=-l gmock      #  libgmock-dev
#LIBS+=-l benchmark  #  libbenchmark-dev

rb_test: $(TESTS) $(SOURCES) $(HEADERS)
#  libgmock-dev libbenchmark-dev
	@echo "Install Google test library: sudo apt-get install libgtest-dev"
	@-rm rb_test 2> /dev/null
	g++ --std=c++20 -Wall -Wextra -Werror -pedantic -O3 -flto=auto -DNDEBUG $(INCLUDE_DIRS) $(TESTS) $(SOURCES) $(LIBS) -o rb_test

test: rb_test
	./rb_test

rb_coverage:  $(TESTS) $(SOURCES) $(HEADERS) coverage-tool
	@echo "Install Google test library: sudo apt-get install libgtest-dev"
	@-rm *.gcno *.gcda 2> /dev/null
	g++ --std=c++20 -Wall -Wextra -Werror -pedantic -Og -g -finline-limit=0 -DNDEBUG $(INCLUDE_DIRS) --coverage -fprofile-abs-path $(TESTS) $(SOURCES) $(LIBS) -o rb_coverage

rb_coverage.info: rb_coverage coverage-tool
	@echo "Install lcov: sudo apt-get install lcov"
	@-rm *.gcda rb_coverage.info 2> /dev/null
	./rb_coverage
	lcov --base-directory . --no-external --rc lcov_branch_coverage=1 --directory . --exclude "*.cc" --capture --output-file rb_coverage.info

rb_cleaned.info: rb_coverage.info coverage-tool
	@-rm rb.info 2> /dev/null
	lcov --rc lcov_branch_coverage=1 --add-tracefile rb_coverage.info --exclude "*.cc" --output-file rb_cleaned.info

coverage/index.html: rb_cleaned.info coverage-tool
	@-rm -rf coverage/ 2> /dev/null
	genhtml --show-details --demangle-cpp --branch-coverage --prefix $(realpath ..) rb_cleaned.info --output-directory coverage/

coverage: rb_coverage rb_cleaned.info
	@echo ""
	@echo "See detailed coverage report at file " $(realpath ./coverage/index.html)

doc/html/index.html: Doxyfile $(SOURCES) $(HEADERS)
	doxygen Doxyfile

doc: doxygen-tool doc/html/index.html $(SOURCES) $(HEADERS)
	@echo "See detailed documentation at file " $(realpath ./doc/html/index.html)

all: test coverage doc

clean:
	@-rm -rf rb_test rb_coverage *.gcno *.gcda *.info coverage/ doc/ 2> /dev/null

format: format-tool $(TESTS) $(SOURCES) $(HEADERS)
	@echo "Install clang-format: sudo apt-get install clang-format"
	clang-format -i $(SOURCES) $(HEADERS) $(TESTS)

cppcheck: cppcheck-tool $(TESTS) $(SOURCES) $(HEADERS)
	cppcheck --platform=unix64 --language=c++ --std=c++20 $(SOURCES) $(HEADERS)
	cppcheck --platform=unix64 --language=c++ --std=c++20 $(TESTS)

.PHONY: format test coverage doc clean


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

.PHONY: cppcheck-tool coverage-tool doxygen-tool format-tool
