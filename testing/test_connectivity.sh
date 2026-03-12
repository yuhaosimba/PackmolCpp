#!/bin/bash
echo "Running connectivity test..."
packmol_bin="${PACKMOL_BIN:-../packmol}"

# Test 1
"$packmol_bin" < ./input_files/benzene2.inp > /dev/null
if diff <(grep "^CONECT" ./output.pdb) <(grep "^CONECT" ./output_files/benzene2.pdb) >/dev/null; then
    exit 0  # Files are identical
else
    echo "Connectivity test failed (test 1)."
    exit 1  # Files are different
fi

# Test 2
"$packmol_bin" < ./input_files/water_box_conect.inp > /dev/null
if diff <(grep "^CONECT" ./output.pdb) <(grep "^CONECT" ./output_files/water_connect.pdb) >/dev/null; then
    exit 0  # Files are identical
else
    echo "Connectivity test failed (test 2)."
    exit 1  # Files are different
fi
