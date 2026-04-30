#!/bin/bash
# Rewrites bare-basename dylib deps in $1 to absolute paths and re-signs.
# Needed because some deal.II-candi / OpenCASCADE libs have bare-basename
# install names, which dyld cannot find on macOS without DYLD_LIBRARY_PATH.
set -e
target="$1"
if [ -z "$target" ] || [ ! -f "$target" ]; then
  echo "fix_dylib_deps.sh: target '$target' not found" >&2
  exit 1
fi

search_dirs=(
  /usr/local/lib
  /Users/saeed/dealii-candi/parmetis-4.0.3/lib
)

for lib in $(otool -L "$target" | awk '/^\t/ {print $1}' | grep -vE '^(/|@)'); do
  for dir in "${search_dirs[@]}"; do
    if [ -f "$dir/$lib" ]; then
      install_name_tool -change "$lib" "$dir/$lib" "$target"
      echo "  patched: $lib -> $dir/$lib"
      break
    fi
  done
done

codesign --force --sign - "$target"
