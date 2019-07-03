astyle  --style=allman  --style=break --indent=tab --indent=force-tab=8 --indent-classes --indent-switches --indent-cases --indent-namespaces --indent-after-parens --indent-preproc-define --indent-preproc-cond --indent-col1-comments --break-blocks --pad-oper --pad-comma --pad-header  --align-pointer=type  --align-reference=type --break-one-line-headers --add-braces  --attach-return-type --remove-comment-prefix --mode=c --suffix=none --preserve-date --verbose --lineend=linux --max-code-length=100  $(pwd)/*.cpp,*.h,*.hpp

cppcheck --enable=all eigen_bicgstab.cpp --report-progress
cppcheck --enable=all --language=c++ bicgstabell.h --report-progress
rm bicgstabell_test.out
g++ eigen_bicgstab.cpp -std=c++1z -I/usr/include/eigen3 -I$(pwd)/eigen-eigen-323c052e1731 -g -Wall -O3 -D PRECONDITIONER="Eigen::DiagonalPreconditioner<T>" -D ALGORITHM=2 -D DEBUG=0 -D ELL=4 -D MAJORITY=Eigen::RowMajor -o bicgstabell_test.out
#g++ eigen_bicgstab.cpp -std=c++1z -I$(pwd)/eigen-eigen-323c052e1731 -g -Wall -O3 -D PRECONDITIONER="Eigen::DiagonalPreconditioner<T>" -D ALGORITHM=2 -D DEBUG=0 -D ELL=4 -D MAJORITY=Eigen::RowMajor -o bicgstabell_test

echo "Done compiling"
./bicgstabell_test.out debug_matrix/input_A.txt debug_matrix/input_b.txt debug_matrix/output.txt 0 12
