#!/bin/bash 

if [[ "$PWD" == */scripts ]]; then
	echo "Script must be run from the project's root directory"
	exit 1
fi

for directory in "kagen" "app"; do
	find $directory -type f \( -name "*.cpp" -o -name "*.cc" -o -name "*.c" -o -name "*.h" \) -exec clang-format -i {} \;
done
