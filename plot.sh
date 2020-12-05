RELATIVE_DIR=$(dirname "$0")

cd $RELATIVE_DIR/tsp_view
cat ../build/result.tour |  python tsp_view.py
cd -