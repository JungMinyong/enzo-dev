#!/bin/bash
set -euo pipefail

# Non-recursive safer macro fixer:
#  - Replace max( -> MAX_VAL( and min( -> MIN_VAL(
#    but DO NOT touch std::max / std::min.
#  - Fix %"ISYM" / %"PSYM" / %"ESYM" / %"FSYM"
#  - Only modifies files in the current directory
#  - Creates .bak backups

echo "Fixing macros ONLY in current directory..."

for file in *.c *.C *.cc *.cpp *.h; do
    # Skip patterns when no files match
    [[ ! -e "$file" ]] && continue

    echo "Processing: $file"

    sed -i.bak \
        -e 's/std::[[:space:]]*max[[:space:]]*(/STD__MAX__(/g' \
        -e 's/std::[[:space:]]*min[[:space:]]*(/STD__MIN__(/g' \
        -e 's/\<max(/MAX_VAL(/g' \
        -e 's/\<min(/MIN_VAL(/g' \
        -e 's/STD__MAX__(/std::max(/g' \
        -e 's/STD__MIN__(/std::min(/g' \
        -e 's/%"\(ISYM\)/%" \1/g' \
        -e 's/%"\(PSYM\)/%" \1/g' \
        -e 's/%"\(ESYM\)/%" \1/g' \
        -e 's/%"\(FSYM\)/%" \1/g' \
        -e 's/%"\(PISYM\)/%" \1/g' \
        -e 's/%"\(GSYM\)/%" \1/g' \
        -e 's/%"\(GOUTSYM\)/%" \1/g' \
        "$file"
done

echo "Done. Backups saved as *.bak."
