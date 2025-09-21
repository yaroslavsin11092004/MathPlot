function build_project()
{
	cmake -DCMAKE_CXX_COMPILER=/usr/bin/clang++ .. 
	make -j${nproc}
	find -type f ! -name "*.so" -delete
	rm -rf CMakeFiles
	rm -rf .cmake
}
if [ -d "./build" ]; then
	cd ./build 
	build_project
else
	mkdir build
	cd build
	build_project
fi
