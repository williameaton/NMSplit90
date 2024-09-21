# Check for flags
skip_cleanup=false
skip_pytest=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --nc) skip_cleanup=true ;;
        --np) skip_pytest=true ;;
    esac
    shift
done



echo "Running test_ylm_values..."
./test_ylm_values
echo "Running test_ylm_function..."
./test_ylm_function


if [ "$skip_pytest" = false ]; then
    pytest ylm/test_ylm.py
fi
if [ "$skip_cleanup" = false ]; then
    cd ylm
    bash cleanup.sh
fi