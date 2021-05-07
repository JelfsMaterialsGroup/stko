STK_DIR=$(readlink -e "$0")
STKO_DIR=${STKO_DIR%/*/*}
rm -- "$STKO_DIR"/docs/source/{stko.*.rst,stk.rst}
sphinx-apidoc -feEM -o "$STKO_DIR/docs/source" "$STKO_DIR/src/stko"