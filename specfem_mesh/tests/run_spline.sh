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


echo "Running test_spline..."
./test_spline


if [ "$skip_pytest" = false ]; then
    pytest spline/test_spline.py
fi


if [ "$skip_cleanup" = false ]; then
    cd spline
    bash cleanup.sh
fi
