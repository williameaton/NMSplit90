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

echo "Running test_mode_norm..."
./test_mode_norm


if [ "$skip_pytest" = false ]; then
    pytest mode_normalisation/test_norm.py
fi
if [ "$skip_cleanup" = false ]; then
    cd mode_normalisation/
    bash cleanup.sh
fi