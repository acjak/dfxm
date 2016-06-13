# find output/ -name '*.png' -o -type d -printf '%h\n' | sort | uniq > non-empty-dirs.tmp
find output/ -name '*.png' -or -name '*.pdf' -or -name "*.npy" -exec dirname {} \; | sort -u > non-empty-dirs.tmp

find output/ -type d -print | sort | uniq > all-dirs.tmp
comm -23 all-dirs.tmp non-empty-dirs.tmp > dirs-to-be-deleted.tmp

less dirs-to-be-deleted.tmp

# cat dirs-to-be-deleted.tmp | xargs rm -rf

find output/ -name '*.inf' -and \! -name "*.inf" -exec dirname {} \; | sort -u
