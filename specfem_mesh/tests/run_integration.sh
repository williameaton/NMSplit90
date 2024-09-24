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

echo "Running test_1_integration..."
./test_1_integration
echo "Running test_ylm_integration..."
./test_ylm_integration


if [ "$skip_pytest" = false ]; then
    pytest integration/test_integration.py
fi
if [ "$skip_cleanup" = false ]; then
    cd integration  
    bash cleanup.sh
fi