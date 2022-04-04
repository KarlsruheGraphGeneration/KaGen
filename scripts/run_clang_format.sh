#!/bin/sh 
if [[ "$PWD" == */scripts ]]; then
	echo "Script must be run from the project's root directory"
	exit 1
fi

for directory in "include" "library" "app"; do
	find $directory -type f \( -name "*.cpp" -o -name "*.cc" -o -name "*.c" -o -name "*.h" \) -exec clang-format -i {} \;
done
