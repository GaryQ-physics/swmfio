print('\n\nimporting named_var_indexes\n\n')

_x                         = 0
_y                         = 1
_z                         = 2
_rho                       = 3
_ux                        = 4
_uy                        = 5
_uz                        = 6
_e                         = 7
_bx                        = 8
_by                        = 9
_bz                        = 10
_b1x                       = 11
_b1y                       = 12
_b1z                       = 13
_p                         = 14
_jx                        = 15
_jy                        = 16
_jz                        = 17 ; nVarNeeded=18;
_norm_b                    = 18
_norm_b1                   = 19
_norm_j                    = 20
_div_b                     = 21
_div_b1                    = 22
_div_j                     = 23
_curl_b_x                  = 24
_curl_b_y                  = 25
_curl_b_z                  = 26
_norm_curl_b               = 27
_curl_b1_x                 = 28
_curl_b1_y                 = 29
_curl_b1_z                 = 30
_norm_curl_b1              = 31
_curl_j_x                  = 32
_curl_j_y                  = 33
_curl_j_z                  = 34
_norm_curl_j               = 35
_relative_div_b            = 36
_relative_div_b1           = 37
_relative_div_j            = 38
_FNdel_b                   = 39
_FNdel_b1                  = 40
_FNdel_j                   = 41
_normalized_div_b          = 42
_normalized_div_b1         = 43
_normalized_div_j          = 44
_div_over_curl_b           = 45
_div_over_curl_b1          = 46
_div_over_curl_j           = 47
_jRx                       = 48
_jRy                       = 49
_jRz                       = 50
_norm_jR                   = 51
_jR_error                  = 52
_jR_fractional_error       = 53
_gridspacing               = 54 ; nVarTot = 55;

index2str = {
    _x                   : 'x'                   ,
    _y                   : 'y'                   ,
    _z                   : 'z'                   ,
    _rho                 : 'rho'                 ,
    _ux                  : 'ux'                  ,
    _uy                  : 'uy'                  ,
    _uz                  : 'uz'                  ,
    _e                   : 'e'                   ,
    _bx                  : 'bx'                  ,
    _by                  : 'by'                  ,
    _bz                  : 'bz'                  ,
    _b1x                 : 'b1x'                 ,
    _b1y                 : 'b1y'                 ,
    _b1z                 : 'b1z'                 ,
    _p                   : 'p'                   ,
    _jx                  : 'jx'                  ,
    _jy                  : 'jy'                  ,
    _jz                  : 'jz'                  ,
    _norm_b              : 'norm_b'              ,
    _norm_b1             : 'norm_b1'             ,
    _norm_j              : 'norm_j'              ,
    _div_b               : 'div_b'               ,
    _div_b1              : 'div_b1'              ,
    _div_j               : 'div_j'               ,
    _curl_b_x            : 'curl_b_x'            ,
    _curl_b_y            : 'curl_b_y'            ,
    _curl_b_z            : 'curl_b_z'            ,
    _norm_curl_b         : 'norm_curl_b'         ,
    _curl_b1_x           : 'curl_b1_x'           ,
    _curl_b1_y           : 'curl_b1_y'           ,
    _curl_b1_z           : 'curl_b1_z'           ,
    _norm_curl_b1        : 'norm_curl_b1'        ,
    _curl_j_x            : 'curl_j_x'            ,
    _curl_j_y            : 'curl_j_y'            ,
    _curl_j_z            : 'curl_j_z'            ,
    _norm_curl_j         : 'norm_curl_j'         ,
    _relative_div_b      : 'relative_div_b'      ,
    _relative_div_b1     : 'relative_div_b1'     ,
    _relative_div_j      : 'relative_div_j'      ,
    _FNdel_b             : 'FNdel_b'             ,
    _FNdel_b1            : 'FNdel_b1'            ,
    _FNdel_j             : 'FNdel_j'             ,
    _normalized_div_b    : 'normalized_div_b'    ,
    _normalized_div_b1   : 'normalized_div_b1'   ,
    _normalized_div_j    : 'normalized_div_j'    ,
    _div_over_curl_b     : 'div_over_curl_b'     ,
    _div_over_curl_b1    : 'div_over_curl_b1'    ,
    _div_over_curl_j     : 'div_over_curl_j'     ,
    _jRx                 : 'jRx'                 ,
    _jRy                 : 'jRy'                 ,
    _jRz                 : 'jRz'                 ,
    _norm_jR             : 'norm_jR'             ,
    _jR_error            : 'jR_error'            ,
    _jR_fractional_error : 'jR_fractional_error' ,
    _gridspacing         : 'gridspacing'         ,
}


str2index = {
    'x'                  :  _x                   ,
    'y'                  :  _y                   ,
    'z'                  :  _z                   ,
    'rho'                :  _rho                 ,
    'ux'                 :  _ux                  ,
    'uy'                 :  _uy                  ,
    'uz'                 :  _uz                  ,
    'e'                  :  _e                   ,
    'bx'                 :  _bx                  ,
    'by'                 :  _by                  ,
    'bz'                 :  _bz                  ,
    'b1x'                :  _b1x                 ,
    'b1y'                :  _b1y                 ,
    'b1z'                :  _b1z                 ,
    'p'                  :  _p                   ,
    'jx'                 :  _jx                  ,
    'jy'                 :  _jy                  ,
    'jz'                 :  _jz                  ,
    'norm_b'             :  _norm_b              ,
    'norm_b1'            :  _norm_b1             ,
    'norm_j'             :  _norm_j              ,
    'div_b'              :  _div_b               ,
    'div_b1'             :  _div_b1              ,
    'div_j'              :  _div_j               ,
    'curl_b_x'           :  _curl_b_x            ,
    'curl_b_y'           :  _curl_b_y            ,
    'curl_b_z'           :  _curl_b_z            ,
    'norm_curl_b'        :  _norm_curl_b         ,
    'curl_b1_x'          :  _curl_b1_x           ,
    'curl_b1_y'          :  _curl_b1_y           ,
    'curl_b1_z'          :  _curl_b1_z           ,
    'norm_curl_b1'       :  _norm_curl_b1        ,
    'curl_j_x'           :  _curl_j_x            ,
    'curl_j_y'           :  _curl_j_y            ,
    'curl_j_z'           :  _curl_j_z            ,
    'norm_curl_j'        :  _norm_curl_j         ,
    'relative_div_b'     :  _relative_div_b      ,
    'relative_div_b1'    :  _relative_div_b1     ,
    'relative_div_j'     :  _relative_div_j      ,
    'FNdel_b'            :  _FNdel_b             ,
    'FNdel_b1'           :  _FNdel_b1            ,
    'FNdel_j'            :  _FNdel_j             ,
    'normalized_div_b'   :  _normalized_div_b    ,
    'normalized_div_b1'  :  _normalized_div_b1   ,
    'normalized_div_j'   :  _normalized_div_j    ,
    'div_over_curl_b'    :  _div_over_curl_b     ,
    'div_over_curl_b1'   :  _div_over_curl_b1    ,
    'div_over_curl_j'    :  _div_over_curl_j     ,
    'jRx'                :  _jRx                 ,
    'jRy'                :  _jRy                 ,
    'jRz'                :  _jRz                 ,
    'norm_jR'            :  _norm_jR             ,
    'jR_error'           :  _jR_error            ,
    'jR_fractional_error':  _jR_fractional_error ,
    'gridspacing'        :  _gridspacing         ,
}
