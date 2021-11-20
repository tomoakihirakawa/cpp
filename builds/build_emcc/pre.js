//Module['TOTAL_MEMORY'] = 65536;
Module['TOTAL_STACK']  = 512 * 1024;

Module['preInit'] = function() {
    ASMJS_PAGE_SIZE  = 65536;
    MIN_TOTAL_MEMORY = 65536;
};
