#ifndef __AC_FIXED_H
#define __AC_FIXED_H
//Catapult-2010a.54beta release
#include <ac_int.h>
//Stuff originally in native ac_int.h but excluded in SLEC version
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> class ac_fixed;
#include <ac_sc.h>
#include <calypto_float_includes.h>
#include <iostream>
#ifdef __ENABLE_ACFIXED_DEBUG_SUPPORT
#include <calypto_synth_system_headers.h>
#endif
#ifndef AC_MAX
#  define AC_MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef AC_MIN
#  define AC_MIN(a,b) ((a) < (b) ? (a) : (b))
#endif


namespace ac_private {
  template<typename T>
  struct rt_ac_fixed_T {
    template<int W, int I, bool S>
    struct op1 {
      typedef typename T::template rt_T< ac_fixed<W,I,S,AC_TRN,AC_WRAP> >::mult mult;
      typedef typename T::template rt_T< ac_fixed<W,I,S,AC_TRN,AC_WRAP> >::plus plus;
      typedef typename T::template rt_T< ac_fixed<W,I,S,AC_TRN,AC_WRAP> >::minus2 minus;
      typedef typename T::template rt_T< ac_fixed<W,I,S,AC_TRN,AC_WRAP> >::minus minus2;
      typedef typename T::template rt_T< ac_fixed<W,I,S,AC_TRN,AC_WRAP> >::logic logic;
      typedef typename T::template rt_T< ac_fixed<W,I,S,AC_TRN,AC_WRAP> >::div2 div;
      typedef typename T::template rt_T< ac_fixed<W,I,S,AC_TRN,AC_WRAP> >::div div2;
    };
  };
}

#pragma calypto_flag DOUBLETOSCBIGINT
sc_bigint<128> double_to_scbigint(double a)
{
       return 0;
}

#pragma calypto_flag CALLREPLACEPRAGMA
template<int W>
bool equal_ones_from(int start, sc_bigint<W> val)
{
    return (~(val >>(start)) == 0) ;
}

template<int W>
bool equal_zeros_to(int start, sc_bigint<W> val)
{
    sc_bigint<W> temp = ~0;
    temp = ~(temp << start);
    return ((val & temp) == 0);
}

template<int N>
inline calypto_double const& ldexpr32(calypto_double& d) {
    calypto_double d2 = d;
    if(N < 0)
    {
        for(int i=0; i < -N; i++)
        {
            d2 /= (calypto_double)((Slong) 1 << 30);
            d2/= 4;
        }
    }
    else
    {
        for(int i=0; i < N; i++)
        {
            d2 = d2 * (calypto_double)((Ulong) 1 << 32);
        }
    }
    d = d2;
    return d;
}

template<> inline calypto_double const& ldexpr32<0>(calypto_double& d) { return d; }

template<> inline calypto_double const& ldexpr32<1>(calypto_double& d)
{
    calypto_double dtmp;
    dtmp = (d * ((calypto_double)((Ulong) 1 << 32)));
    d = dtmp;
    return d;
}

template<> inline calypto_double const& ldexpr32<-1>(calypto_double& d) {
    calypto_double dtmp = d / (calypto_double)((long long) 1 << 30);
    dtmp = dtmp/ (calypto_double)(1 << 2);
    d = dtmp;
    return d;
}

template<> inline calypto_double const& ldexpr32<2>(calypto_double& d)
{
    calypto_double dtmp =
        (d * (calypto_double)((Ulong) 1 << 32));
    dtmp= dtmp * (calypto_double)((Ulong) 1 << 32);
    d = dtmp;
    return d;
}

template<> inline calypto_double const& ldexpr32<-2>(calypto_double& d)
{
    calypto_double dtmp =
        (d / (calypto_double)((Ulong) 1 << 30));

    dtmp /=(1 <<2 );
    dtmp /= (calypto_double) ((Ulong) 1 << 30);
    dtmp /=(1 << 2);
    d = dtmp;
    return d;
}

template<int N>
inline calypto_double const& ldexpr( calypto_double& d) {
    if(N < 0 )
    {
        int temp = -N &31;
        int t;
        while (temp > 30 )
        {
            t = (unsigned) 1 << 30;
            temp -= 30 ;
            d = d/t;
        }
        t = (unsigned) 1 << temp;
        d = d/t;
    }
    else
    {
        calypto_double divisor =
            (calypto_double)( (unsigned) 1 << (N & 31));

        d= ( d *divisor);

    }
    ldexpr32<N/32>( d);
    return d;
}

template<int N, int W, bool S>
inline void conv_from_fraction(calypto_double& d, sc_bigint<W> *r, bool *qb, bool *rbits, bool *o) {
    bool b = d < 0;
    calypto_double d2 = b ? -d : d;
    // calypto_double dfloor = mgc_floor(d2); NIKHIL
//    calypto_double dfloor = 0;
    *o = false; // NIKHIL
    //d2 = d2 ;//- dfloor;
    *r = 0;
    for(int i=N-1; i >=0; i--) {
      *r <<= 32;
      d2 = d2 * (calypto_double) ((Ulong) 1 << 32);
      unsigned long long k_temp =
          float64_to_uint64_round_to_zero(d2.data());
      unsigned k = static_cast<unsigned int>( k_temp);
      *r |= b ? ~k : k;
      d2 -= (calypto_double) k;
    }
    d2 = d2 * 2.0;
    signed int k_temp2 = float64_to_int32_round_to_zero(d2.data());
    bool k = (static_cast<int>(k_temp2)) != 0;  // is 0 or 1
    d2 -= (calypto_double)k;
    *rbits = d2 != 0.0;
    *qb = (b && *rbits) ^ k;
    if(b && !*rbits && !*qb)
        *r = *r + 1;
    *o |= b ^ (*r < 0);
  }

  template<int W>
  inline calypto_double to_double_from_sc_bigint(sc_bigint<W> v)
  {
    int N = (W + 31 )/32;
    int shift_ammount = W - W%32;
    sc_bigint<W> temp = v >> shift_ammount;
    calypto_double a = temp.to_int();
    //int x = a; using internal softfloat function to avoid CPT-STS
    int x = static_cast<int>(float64_to_int32_round_to_zero(a.data_));
    v = v << (W%32);
    for(int i=N-2; i >= 0; i--)
    {
        shift_ammount = W - 32*i;
        temp = v >> shift_ammount;
        a = a *(double)((Slong)1 << 32);
        a += (calypto_double)((unsigned )(temp.to_int()));
        v = v << 32;
    }
    return a;
  }

// In general use sc_bigint<W> for the internal storage
template <int W, bool S>
struct ac_fixed_dt {
    typedef sc_bigint<W> dt;
};

// ... but if we are unsigned, use sc_biguint<W>
template <int W>
struct ac_fixed_dt<W, false> {
    typedef sc_biguint<W> dt;
};

//end ac_int.h
#pragma calypto_flag AC_FIXED_DATA
#pragma calypto_flag DATATYPEPRAGMA
#pragma calypto_flag DONT_FLATTEN
template<int W, int I, bool S=true, ac_q_mode Q=AC_TRN, ac_o_mode O=AC_WRAP>
class ac_fixed
{
    public:
#pragma calypto_flag AC_FIXED_DATA
#pragma calypto_flag DONT_FLATTEN
#pragma calypto_flag DATATYPEPRAGMA
        ac_fixed_dt<W,S>::dt val_;

  static const int width = W;
  static const int i_width = I;
  static const bool sign = S;
  static const ac_o_mode o_mode = O;
  static const ac_q_mode q_mode = Q;

#pragma calypto_flag AC_FIXED_SETVAL
template<ac_special_val V> inline ac_fixed &set_val()
{
    if(V == AC_VAL_DC) {
      sc_lv <W> temp;
      this->val_ = temp;
    }
    else if(V == AC_VAL_0)
    {
        this->val_ = 0;
    }
    else if (V == AC_VAL_QUANTUM)
    {
        this->val_ = 1;
    }
    else if (V == AC_VAL_MIN)
    {
        this->val_ = 0;
        this->val_[W-1] = S;
        if (O==AC_SAT_SYM)
            this->val_[0] = 1;
    }
    else if (V == AC_VAL_MAX)
    {
        this->val_ = ~(ac_fixed_dt<W,S>::dt(0));
        this->val_[W-1] = !S;
    }
    return *this;
}

    template<int W_temp>
    inline bool quantization_adjust(bool qb, bool r, bool s, sc_bigint<W_temp> &temp)
    {
/*       if quantization mode is AC_TRN the function is not called.
         So commenting out this line to avoid useless warnings
         if(Q==AC_TRN)
            return false;*/
        if(Q==AC_RND_ZERO)
            qb &= s || r;
        else if(Q==AC_RND_MIN_INF)
            qb &= r;
        else if(Q==AC_RND_INF)
            qb &= !s || r;
        else if(Q==AC_RND_CONV)
            qb &= (((temp & 1) != 0)) || r;
        else if(Q==AC_RND_CONV_ODD)
            qb &= (((temp & 1) == 0)) || r;
        else if(Q==AC_TRN_ZERO)
            qb = s && ( qb || r );
        sc_bigint<W_temp > old_temp_val = temp;
        sc_biguint<W_temp> temp_ext = temp;
        sc_biguint<W_temp + 1 > temp_ext1 = temp_ext;
        temp_ext1 += qb;
        temp = temp_ext1;
        sc_bigint <1> t = temp_ext1 >> (W_temp );
        if(temp == old_temp_val)
            t = 0;
        return (t != 0);
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    inline void overflow_adjust(bool underflow, bool overflow) {
        // AC_SAT_ZERO overflow
        if (O == AC_SAT_ZERO) {
            // Flip the order of the data branches to mimic catapult
            bool needAdjust = overflow || underflow;
            if (!needAdjust) {
            } else {
                val_ = 0;
            }
        } else {
            // mcase rewrite, compute msb, mid, lsb separately like catapult does
            sc_biguint<1> msb, lsb;
            sc_biguint<W> mid;
            sc_biguint<W> onesString = 0;
            if (W > 2)
                onesString = (sc_biguint<W>(1) << (W-2)) - 1; // a string of W-2 ones.

            // compute the msb
            if (overflow) {
                switch (O) {
                // AC_SAT_SYM and AC_SAT: clamp to the max value
                // For signed this is like   {0, 1111, 1}
                // For unsigned this is like {1, 1111, 1}
                case AC_SAT_SYM:
                case AC_SAT:
                    msb = S ? 0 : 1;
                    break;
                default: // no-op
                    msb = val_[W-1];
                    break;
                }
            } else if (underflow) {
                switch (O) {
                // AC_SAT_SYM: clamp to the min value
                // For signed this is like   {1, 0000, 1}
                // For unsigned this is like {0, 0000, 0}
                case AC_SAT_SYM:
                    msb = S ? 1 : 0;
                    break;

                // AC_SAT: clamp to the min value
                // For signed this is like   {1, 0000, 0}
                // For unsigned this is like {0, 0000, 0}
                case AC_SAT: // clamp to the min value
                    msb = S ? 1 : 0;
                    break;
                default: // no-op
                    msb = val_[W-1];
                    break;
                }
            } else {
                msb = val_[W-1];
            }

            // compute the (mid, lsb).  Note that the branch order is different to mimic catapult
            if (underflow) {
                switch (O) {
                // AC_SAT_SYM: clamp to the min value
                // For signed this is like   {1, 0000, 1}
                // For unsigned this is like {0, 0000, 0}
                case AC_SAT_SYM:
                    mid = 0;
                    lsb = S ? 1 : 0;
                    break;

                // AC_SAT: clamp to the min value
                // For signed this is like   {1, 0000, 0}
                // For unsigned this is like {0, 0000, 0}
                case AC_SAT:
                    mid = 0;
                    lsb = 0;
                    break;
                default: // no-op
                    mid = (val_ >> 1) & onesString;
                    lsb = val_[0];
                    break;
                }
            } else if (overflow) {
                switch (O) {
                // AC_SAT_SYM and AC_SAT: clamp to the max value
                // For signed this is like   {0, 1111, 1}
                // For unsigned this is like {1, 1111, 1}
                case AC_SAT_SYM:
                case AC_SAT:
                    mid = onesString;
                    lsb = 1;
                    break;
                default: // no-op
                    mid = (val_ >> 1) & onesString;
                    lsb = val_[0];
                    break;
                }
            } else {
                mid = (val_ >> 1) & onesString;
                lsb = val_[0];
            }

            if (W == 1) {
                val_ = msb;
            } else if (W == 2) {
                val_ = (msb << 1) | lsb;
            } else {
                val_ = (msb << W-1) | (mid << 1) | lsb;
            }
        }
    }
    inline bool is_neg() const { return S && val_<0;}

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed(ac_fixed_dt<W,S>::dt val_arg) : val_(val_arg) {}

    #pragma calypto_flag AC_FIXED_DEFAULT_CTOR
    ac_fixed() { }


    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed(ac_fixed const& op)
    {
        // invoke operator= so that op can be rounded if necessary
        (*this) = op;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    inline ac_fixed (ac_fixed<W2,I2,S2,Q2,O2> const &op)
    {
        #pragma calypto_flag IGNORE_INIT
        enum {F=W-I, F2=W2-I2,
              // W_temp : wide enough so that "op.val_ << (F-F2)" won't overflow
              // also op.val_ >> (F2 - F) shouldn't overflow in the case that F2-F<0
              W_temp = W2+!S2 + AC_MAX(F-F2,0) + AC_MAX(-1*(F2-F),0),
              QUAN_INC = F2>F && !(Q==AC_TRN || (Q==AC_TRN_ZERO && !S2)) };
        sc_bigint<W_temp> temp ;

        if(F == F2)
            temp = op.val_;
        else if (F2 > F)
        {
            temp = op.val_ >> (F2 - F);
            if(Q!=AC_TRN && !(Q==AC_TRN_ZERO && !S2))
            {

                #pragma calypto_flag IGNORE_INIT
                bool qb;
                if(F2-F > W2)
                    qb =  (op.val_ < 0) ;
                else
                    qb = (bool) op[F2 - F - 1];

                #pragma calypto_flag IGNORE_INIT
                bool r ;
                if(F2 > F+1)
                    r = !(equal_zeros_to<W2>(F2-F-1, op.val_));
                else
                    r = false;
                if(Q!=AC_TRN)
                    quantization_adjust(qb, r, S2 && op.val_ < 0, temp);
            }
        }
        else
        {
            temp = op.val_ << (F-F2);
        }
        val_ = temp;
        // handle overflow/underflow
        if(O!=AC_WRAP && ((!S && S2) || I-S < I2-S2+(QUAN_INC || (S2 && O==AC_SAT_SYM && (O2 != AC_SAT_SYM || F2 > F) ))) ) { // saturation
            // From the System-C spec:
            // The following flags and symbols are used in the template above and in Table 37:
            //   x represents a binary digit (0 or 1).
            //   sD represents a sign bit before overflow handling.
            //   D represents deleted bits.
            //   lD represents the least significant deleted bit.
            //   sR represents the bit on the MSB position of the result number.
            //
            // Overflow shall occur when the value of at least one of the deleted bits (sD, D, lD) is not equal to the original
            // value of the bit on the MSB position of the result (sRo).

            // Reference code: sc_bigint<W_temp> deletedBits = temp >> W
            // Optimization 1: include msb of val_ in temp.  This simplifies some logic below
            // Optimization 2: don't include the msb of temp in deletedBits.  This duplicates some checks below
            // Optimization 3: inline the deletedBits definition.  Storing it to a temporary risks sign extension, etc
            // Optimized code: temp.range(W_temp-2,W-1)

            // Reference code: bool deletedBits_has_a_one = deletedBits != 0;
            // Reference code: bool deletedBits_has_a_zero = ~deletedBits != 0;
            bool deletedBits_has_a_one_or_neg_val;
            bool deletedBits_has_a_zero_or_val_positive;
            if (S) {
                if (W_temp-2 >= W-1) {
                    // compute an op on {deleted bits, sign of val_} = {temp[W_temp-2:W], temp[W-1]}
                    deletedBits_has_a_one_or_neg_val = temp.range(W_temp-2, W-1).or_reduce();
                    deletedBits_has_a_zero_or_val_positive = ! temp.range(W_temp-2,W-1).and_reduce();
                } else {
                    // there are no deleted bits.  Just pick up the sign of val_
                    deletedBits_has_a_one_or_neg_val = (val_ < 0);
                    deletedBits_has_a_zero_or_val_positive = (val_ >= 0);
                }
            } else {
                if (W_temp-2 >= W) {
                    // compute |{deleted bits, sign of val_} = |{temp[W_temp-2:W], 0}
                    deletedBits_has_a_one_or_neg_val = temp.range(W_temp-2,W).or_reduce();
                } else {
                    // there are no deleted bits.  Just pick up the sign of val_
                    deletedBits_has_a_one_or_neg_val = false;
                }

                // val_ is always positive
                deletedBits_has_a_zero_or_val_positive = true;
            }

            // Reference code: bool overflow = (temp >= 0) && (deletedBits_has_a_one || (val_ < 0));
            bool overflow = (temp >= 0) && deletedBits_has_a_one_or_neg_val;

            // Reference code: bool valIsValidNegativeNumber = (val_ < 0);
            // (This is now built into deletedBits_has_a_one_or_neg_val)
            bool valIsValidNegativeNumber = true;
            if (O==AC_SAT_SYM && S && S2 && (W >= 2)) {  // values like 1000000 (INT_MIN) aren't allowed in AC_SAT_SYM mode
                // Reference code: sc_biguint<W> onesString = (sc_biguint<W>(1) << (W-1)) - 1; // a string of W-1 ones.
                // Reference code: valIsValidNegativeNumber &= !((val_ & onesString) == 0);
                valIsValidNegativeNumber &= !(val_.range(W-2,0) == 0);
            }
            // Reference code:  bool underflow = (temp < 0) && (deletedBits_has_a_zero || !valIsValidNegativeNumber || !S);
            bool underflow = (temp < 0) && (deletedBits_has_a_zero_or_val_positive || !valIsValidNegativeNumber);

            overflow_adjust(underflow, overflow);
        }
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =(ac_fixed const& op)
    {
        val_ = op.val_;

        // AC_SAT_SYM disllows the value 1000000.  round to 1000001
//        sc_biguint<W> onesString = (sc_biguint<W>(1) << (W-1)) - 1; // a string of W-1 ones.
//        if(S && O==AC_SAT_SYM) {
//            if((val_ < 0) && ((val_ & onesString) == 0))
//                val_[0] = 1;
//        }

        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed& operator = (ac_fixed<W2,I2,S2,Q2,O2> const &op)
    {
        ac_fixed temp(op);
        this->val_ = temp.val_;
        return *this;
    }


    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, bool S2>
    inline ac_fixed (ac_int<W2,S2> const &op) {
        #pragma calypto_flag IGNORE_INIT
        ac_fixed<W2, W2, S2> f_op;
        f_op.val_ = to_sc(op);
        *this = f_op;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, bool S2>
    ac_fixed& operator = (ac_int<W2,S2> const &op) {
        ac_fixed temp(op);
        this->val_ = temp.val_;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( bool b ) {
        *this = (ac_int<1,false>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( char b )
    {
        *this = (ac_int<8,true>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( signed char b )
    {
        *this = (ac_int<8,true>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( unsigned char b )
    {
        *this = (ac_int<8,false>) b;

    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( signed short b )
    {
        *this = (ac_int<16,true>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( unsigned short b )
    {
        *this = (ac_int<16,false>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( signed int b )
    {
        *this = (ac_int<32,true>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( unsigned int b )
    {
        *this = (ac_int<32,false>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( signed long b )
    {
        *this = (ac_int<32,true>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( unsigned long b )
    {
        *this = (ac_int<32,false>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( Slong b )
    {
        *this = (ac_int<64,true>) b;
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( Ulong b )
    {
        *this = (ac_int<64,false>) b;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator = ( bool b )
    {
        *this = (ac_int<1,false>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( char b )
    {
        *this = (ac_int<8,true>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( signed char b )
    {
        *this = (ac_int<8,true>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( unsigned char b )
    {
        *this = (ac_int<8,false>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( signed short b )
    {
        *this = (ac_int<16,true>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( unsigned short b )
    {
        *this = (ac_int<16,false>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( signed int b )
    {
        *this = (ac_int<32,true>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( unsigned int b )
    {
        *this = (ac_int<32,false>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( signed long b )
    {
        *this = (ac_int<32,true>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( unsigned long b )
    {
        *this = (ac_int<32,false>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( Slong b )
    {
        *this = (ac_int<64,true>) b;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    #pragma calypto_flag CALLREPLACEPRAGMA
    ac_fixed& operator =( Ulong b )
    {
        *this = (ac_int<64,false>) b;
        return *this;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    void double_to_ac_fixed_ctor( const calypto_double& d){

 //       val_ =0;

        calypto_double di =d;
        ldexpr<(-(I+!S+((32-W-!S)&31)))>(di);
        bool o, qb, r;
        bool neg_src = d < 0;
        #pragma calypto_flag IGNORE_INIT
        enum { W_temp = ((W + !S + 31)/32)*32} ;
        sc_bigint<W_temp> temp ;
//        temp = val_;
        conv_from_fraction<(W+31+!S)/32, W_temp, S>(di, &temp, &qb, &r, &o);
        if(Q!=AC_TRN)
            quantization_adjust(qb, r, neg_src, temp);
        val_ = temp;
        // a neg number may become non neg (0) after quantization
        neg_src &= o || temp < 0;

        if(O!=AC_WRAP) { // saturation
            bool overflow, underflow;
            sc_bigint<W+!S> mask = 1;
            mask <<= (W-1);
            //bool neg_trg = S &&((val_ & mask)!=0);
            bool neg_trg = S &&(val_[W-1]!= 0);
            if(o)
            {
                overflow = !neg_src;
                underflow = neg_src;
            }
            else
            {
            //    bool deleted_bits_zero = !((temp >> W)==0);
                int temp_int = (temp>> W).to_int();
                bool deleted_bits_one =   !(~temp_int);
                bool deleted_bits_zero =  !temp_int;

                overflow = !neg_src && (neg_trg || !deleted_bits_zero);
                underflow = neg_src && (!neg_trg || !deleted_bits_one);
            }
            if(O==AC_SAT_SYM && S)
                underflow |= neg_src && (W > 1 ?
                        equal_zeros_to<W>(W-1, val_) : true);

            overflow_adjust(underflow, overflow);
         }
    }

    #pragma calypto_flag AC_FIXED_CTOR
    #pragma calypto_flag CALLREPLACEPRAGMA
    inline ac_fixed( const calypto_double& d ) {
        double_to_ac_fixed_ctor(d);
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag CALYPTO_FLOAT_CONST_OP
    inline ac_fixed( const double d_val ) {
        ac_fixed<128,64> tmp;
        tmp.val_ = double_to_scbigint(d_val);
        ac_fixed tmp1 = ac_fixed ( tmp);
        this->val_ = tmp1.val_;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag CALYPTO_FLOAT_CONST_OP
    ac_fixed operator =( const double d_val ) {
        ac_fixed<128,64> tmp;
        tmp.val_ = double_to_scbigint(d_val);
        ac_fixed tmp1 = ac_fixed ( tmp);
        this->val_ = tmp1.val_;
        return *this;
    }

  // returns false if number is denormal
  template<int WE, bool SE>
  bool normalize(ac_int<WE,SE> &exp) {
    ac_int<W,S> m = this->template slc<W>(0);
    bool r = m.normalize(exp);
    this->set_slc(0,m);
    return r;
  }

    // Explicit conversion functions to ac_int that captures all integer bits (bits are truncated)
    inline ac_int<AC_MAX(I,1),S> to_ac_int() const
    {
        ac_int<AC_MAX(I,1),S> temp = (ac_int<AC_MAX(I,1),S> )(((ac_fixed<AC_MAX(I,1),AC_MAX(I,1),S>) *this).val_);
        return temp;
    }

    inline int to_int() const
    {
        return ((I-W) >= 32) ? 0 : (signed int) (((ac_fixed<AC_MAX(I,1),AC_MAX(I,1),S>) *this).val_.to_int64());
    }

    inline unsigned to_uint() const
    {
        return ((I-W) >= 32) ? 0 : (unsigned int) (((ac_fixed<AC_MAX(I,1),AC_MAX(I,1),S>) *this).val_.to_int64());
    }

    inline long to_long() const
    {
        return ((I-W) >= 32) ? 0 : (signed long) (((ac_fixed<AC_MAX(I,1),AC_MAX(I,1),S>) *this).val_.to_int64());
    }

    inline unsigned long to_ulong() const
    {
        return ((I-W) >= 32) ? 0 : (unsigned long)(((ac_fixed<AC_MAX(I,1),AC_MAX(I,1),S>) *this).val_.to_int64());
    }

    inline Slong to_int64() const
    {
        return ((I-W) >= 64) ? 0  :
            (Slong) (((ac_fixed<AC_MAX(I,1),AC_MAX(I,1),S>) *this).val_.to_int64());
    }

    inline Ulong to_uint64() const
    {
        return ((I-W) >= 64) ? 0 :
            (Ulong) (((ac_fixed<AC_MAX(I,1),AC_MAX(I,1),S>) *this).val_.to_int64());
    }

  #pragma calypto_flag API_NOT_SUPPORTED
    inline double to_double() const
    {
        double d1 = 0.0;
        return d1;
    }

    inline int length() const { return W; }

  template<int W2, int I2, bool S2>
  struct rt {
    enum {
      F=W-I,
      F2=W2-I2,
      mult_w = W+W2,
      mult_i = I+I2,
      mult_s = S||S2,
      plus_w = AC_MAX(I+(S2&&!S),I2+(S&&!S2))+1+AC_MAX(F,F2),
      plus_i = AC_MAX(I+(S2&&!S),I2+(S&&!S2))+1,
      plus_s = S||S2,
      minus_w = AC_MAX(I+(S2&&!S),I2+(S&&!S2))+1+AC_MAX(F,F2),
      minus_i = AC_MAX(I+(S2&&!S),I2+(S&&!S2))+1,
      minus_s = true,
      div_w = W+AC_MAX(W2-I2,0)+S2,
      div_i = I+(W2-I2)+S2,
      div_s = S||S2,
      logic_w = AC_MAX(I+(S2&&!S),I2+(S&&!S2))+AC_MAX(F,F2),
      logic_i = AC_MAX(I+(S2&&!S),I2+(S&&!S2)),
      logic_s = S||S2
    };
    typedef ac_fixed<mult_w, mult_i, mult_s> mult;
    typedef ac_fixed<plus_w, plus_i, plus_s> plus;
    typedef ac_fixed<minus_w, minus_i, minus_s> minus;
    typedef ac_fixed<logic_w, logic_i, logic_s> logic;
    typedef ac_fixed<div_w, div_i, div_s> div;
    typedef ac_fixed<W, I, S> arg1;
  };

  template<typename T>
  struct rt_T {
    typedef typename ac_private::map<T>::t map_T;
    typedef typename ac_private::rt_ac_fixed_T<map_T>::template op1<W,I,S>::mult mult;
    typedef typename ac_private::rt_ac_fixed_T<map_T>::template op1<W,I,S>::plus plus;
    typedef typename ac_private::rt_ac_fixed_T<map_T>::template op1<W,I,S>::minus minus;
    typedef typename ac_private::rt_ac_fixed_T<map_T>::template op1<W,I,S>::minus2 minus2;
    typedef typename ac_private::rt_ac_fixed_T<map_T>::template op1<W,I,S>::logic logic;
    typedef typename ac_private::rt_ac_fixed_T<map_T>::template op1<W,I,S>::div div;
    typedef typename ac_private::rt_ac_fixed_T<map_T>::template op1<W,I,S>::div2 div2;
    typedef ac_fixed<W, I, S> arg1;
  };

  struct rt_unary {
    enum {
      neg_w = W+1,
      neg_i = I+1,
      neg_s = true,
      mag_sqr_w = 2*W-S,
      mag_sqr_i = 2*I-S,
      mag_sqr_s = false,
      mag_w = W+S,
      mag_i = I+S,
      mag_s = false
    };
    typedef ac_fixed<neg_w, neg_i, neg_s> neg;
    typedef ac_fixed<mag_sqr_w, mag_sqr_i, mag_sqr_s> mag_sqr;
    typedef ac_fixed<mag_w, mag_i, mag_s> mag;
    template<unsigned N>
    struct set {
      enum { sum_w = W + ac::log2_ceil<N>::val, sum_i = (sum_w-W) + I, sum_s = S};
      typedef ac_fixed<sum_w, sum_i, sum_s> sum;
    };
  };

// Arithmetic : Binary ----------------------------------------------------
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    typename rt<W2,I2,S2>::mult operator *( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        typename rt<W2,I2,S2>::mult r(this->val_ * op2.val_);
        return r;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    typename  rt<W2,I2,S2>::plus operator +( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };
        #pragma calypto_flag IGNORE_INIT
        typename rt<W2,I2,S2>::plus r;
        if(F == F2)
            r.val_ = this->val_ + op2.val_;
        else if (F > F2)
            r.val_ = (op2.val_ << (F - F2)) + this->val_;
        else
            r.val_ = (this->val_ << (F2 - F)) + op2.val_;
        return r;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    typename rt<W2,I2,S2>::minus operator -( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };
        #pragma calypto_flag IGNORE_INIT
        typename rt<W2,I2,S2>::minus r;
        if(F == F2)
            r.val_ = (this->val_ - op2.val_);
        else if( F > F2)
            r.val_ = this->val_ - (op2.val_ << (F - F2 ) ) ;
        else
            r.val_ = (this->val_ << F2 - F) - op2.val_;
        return r;
    }

    #pragma calypto_flag AC_FIXED_DIVIDE
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    typename rt<W2,I2,S2>::div operator /( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        typename rt<W2,I2,S2>::div r;
        // this is not quite right.  require futher testing
//        if (op2.val_ == 0) {
//            static ac_fixed<W,I,S> max;   // Note: this is max wrt this type, not the return type
//            max.set_val<AC_VAL_MAX>();
//            r = max;
//        } else {
            ac_fixed<W+AC_MAX(W2-I2,0), I, S> t(*this) ;
            r.val_ = t.val_ / op2.val_;
//        }
        return r;
    }
    // Arithmetic assign  ------------------------------------------------------
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed &operator *=( const ac_fixed<W2,I2,S2,Q2,O2> &op2) {
        *this = this->operator *(op2);
        return *this;
    }
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed &operator +=( const ac_fixed<W2,I2,S2,Q2,O2> &op2) {
        *this = this->operator +(op2);
        return *this;
    }
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed &operator -=( const ac_fixed<W2,I2,S2,Q2,O2> &op2)
    {
        *this = this->operator -(op2);
        return *this;
    }
    #pragma calypto_flag AC_FIXED_DIVEQUALS
    #pragma calypto_flag CALLREPLACEPRAGMA
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed &operator /=( const ac_fixed<W2,I2,S2,Q2,O2> &op2)
    {
        *this = this->operator /(op2);
        return *this;
    }
    // increment/decrement by quantum (smallest difference that can be represented)
    // Arithmetic prefix increment, decrement ---------------------------------
    ac_fixed &operator ++()
    {
        ac_fixed<1,I-W+1,false> q;
        q.val_ = 1 ;
        operator += (q);
        return *this;
    }
    ac_fixed &operator --() {
        ac_fixed<1,I-W+1,false> q;
        q.val_ = 1 ;
        operator -= (q);
        return *this;
    }
    // Arithmetic postfix increment, decrement ---------------------------------
    const ac_fixed operator ++(int)
    {
        ac_fixed t = *this;
        ac_fixed<1,I-W+1,false> q;
        q.val_ = 1 ;
        operator += (q);

        return t;
    }
    const ac_fixed operator --(int)
    {
        ac_fixed t = *this;
        //operator-=((ac_fixed<1,I-W+1,false>) 1);
        ac_fixed<1,I-W+1,false> q;
        q.val_ = 1 ;
        operator -= (q);
        return t;
    }
    // Arithmetic Unary --------------------------------------------------------
    ac_fixed operator +()
    {
        return *this;
    }
    ac_fixed<W+1,I+1,true> operator -() const
    {
        return (((ac_fixed<1,1,false>) 0) - *this);
    }
    bool operator ! () const {
        return (val_ == 0 ) ;
    }

    // Bitwise (not arithmetic) unary: complement  -----------------------------
    ac_fixed<W+!S, I+!S, true> operator ~() const
    {
        ac_fixed<W+!S,I+!S,true> r;
        r.val_ = this->val_;
        r.val_ = ~r.val_;
        return r;
    }
    
    ac_fixed<W, I, false> bit_complement() const 
    {
        ac_fixed<W, I, false> r;
        r.val_ = this->val_;
        r.val_ = ~r.val_;
        return r;
    }

    // Bitwise (not arithmetic): and, or, xor ----------------------------------
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    typename rt<W2,I2,S2>::logic operator &( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };

        // Note: it is necessary to cast lhs and rhs up to the correct bit width.
        // Else it is possible for:
        //   ac_fixed<32,32,true> a = 1111111....,   ac_fixed<32,32,false> b = XXXXX
        //   a | b = 1111111
        //   output rt is 33 bits, so cast to 33 bit output:  0111111
        //   this is always > 0, hence underflow will be averted

        typename rt<W2,I2,S2>::logic lhs, rhs;
        if(F == F2)
        {
            lhs.val_ = this->val_;
            rhs.val_ = op2.val_;
        }
        else if(F > F2)
        {
            lhs.val_ = this->val_;
            rhs.val_ = op2.val_ << (F-F2);
        }
        else
        {
            lhs.val_ = this->val_ << (F2 - F);
            rhs.val_ = op2.val_;
        }

        typename rt<W2,I2,S2>::logic r;
        r.val_ = lhs.val_ & rhs.val_;
        return r;
    }

    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    typename rt<W2,I2,S2>::logic operator |( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };

        // Note: it is necessary to cast lhs and rhs up to the correct bit width.
        // Else it is possible for:
        //   ac_fixed<32,32,true> a = 1111111....,   ac_fixed<32,32,false> b = XXXXX
        //   a | b = 1111111
        //   output rt is 33 bits, so cast to 33 bit output:  0111111
        //   this is always > 0, hence underflow will be averted

        typename rt<W2,I2,S2>::logic lhs, rhs;
        if(F == F2)
        {
            lhs.val_ = this->val_;
            rhs.val_ = op2.val_;
        }
        else if(F > F2)
        {
            lhs.val_ = this->val_;
            rhs.val_ = op2.val_ << (F-F2);
        }
        else
        {
            lhs.val_ = this->val_ << (F2 - F);
            rhs.val_ = op2.val_;
        }

        typename rt<W2,I2,S2>::logic r;
        r.val_ = lhs.val_ | rhs.val_;
        return r;
    }
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    typename rt<W2,I2,S2>::logic operator ^( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };

        // Note: it is necessary to cast lhs and rhs up to the correct bit width.
        // Else it is possible for:
        //   ac_fixed<32,32,true> a = 1111111....,   ac_fixed<32,32,false> b = XXXXX
        //   a | b = 1111111
        //   output rt is 33 bits, so cast to 33 bit output:  0111111
        //   this is always > 0, hence underflow will be averted

        typename rt<W2,I2,S2>::logic lhs, rhs;
        if(F == F2)
        {
            lhs.val_ = this->val_;
            rhs.val_ = op2.val_;
        }
        else if(F > F2)
        {
            lhs.val_ = this->val_;
            rhs.val_ = op2.val_ << (F-F2);
        }
        else
        {
            lhs.val_ = this->val_ << (F2 - F);
            rhs.val_ = op2.val_;
        }

        typename rt<W2,I2,S2>::logic r;
        r.val_ = lhs.val_ ^ rhs.val_;
        return r;
    }

    // Bitwise assign (not arithmetic): and, or, xor ----------------------------
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed &operator &= ( const ac_fixed<W2,I2,S2,Q2,O2> &op2 )
    {
        *this = this->operator &(op2);
        return *this;
    }
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed &operator |= ( const ac_fixed<W2,I2,S2,Q2,O2> &op2 )
    {
        *this = this->operator |(op2);
        return *this;
    }
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    ac_fixed &operator ^= ( const ac_fixed<W2,I2,S2,Q2,O2> &op2 )
    {
        *this = this->operator ^(op2);
        return *this;
    }

    // Shift (result constrained by left operand) -------------------------------
    #pragma calypto_flag AC_FIXED_S_LSHIFT
    template<int W2>
    ac_fixed<W,I,S> operator << ( const ac_int<W2,true> &op2 ) const
    {
        enum { ShiftWidth = AC_MIN(W2, 32) };   // Catapult considers at most 32 bits

        ac_fixed r;
        sc_bigint<ShiftWidth> shift = to_sc(op2);
        if(shift < 0 )
        {
            r.val_ = this->val_ >> -shift;
        }
        else
        {
            r.val_ = this->val_ << shift;
        }

        return r;
    }

    #pragma calypto_flag AC_FIXED_LSHIFT
    template<int W2>
    ac_fixed<W,I,S> operator << ( const ac_int<W2,false> &op2 ) const
    {
        enum { ShiftWidth = AC_MIN(W2, 32) };   // Catapult considers at most 32 bits
        sc_biguint<ShiftWidth> shift = to_sc(op2);

        ac_fixed r;
        r.val_ = this->val_ << shift;
        return r;
    }
    #pragma calypto_flag AC_FIXED_S_RSHIFT
    template<int W2>
    ac_fixed<W,I,S> operator >> ( const ac_int<W2,true> &op2 ) const
    {
        enum { ShiftWidth = AC_MIN(W2, 32) };   // Catapult considers at most 32 bits

        ac_fixed r;
        sc_bigint<ShiftWidth> shift = to_sc(op2);
        if(shift < 0 )
        {
            r.val_ = this->val_ << -shift;
        }
        else
        {
            r.val_ = this->val_ >> shift;
        }
        return r;
    }
    #pragma calypto_flag AC_FIXED_RSHIFT
    template<int W2>
    ac_fixed<W,I,S> operator >> ( const ac_int<W2,false> &op2 ) const
    {
        enum { ShiftWidth = AC_MIN(W2, 32) };   // Catapult considers at most 32 bits
        sc_biguint<ShiftWidth> shift = to_sc(op2);

        ac_fixed r;
        r.val_ = this->val_ >> shift;
        return r;
    }

    #pragma calypto_flag AC_FIXED_LSHIFT_ASSIGN
    template<int W2, bool S2>
    ac_fixed &operator <<= ( const ac_int<W2,S2> &op2 ) {
        ac_fixed<W,I,S> r = this->operator<<( op2).val_;
        *this = r;
        return *this;
    }

    #pragma calypto_flag AC_FIXED_RSHIFT_ASSIGN
    template<int W2, bool S2>
    ac_fixed &operator >>= ( const ac_int<W2,S2> &op2 ) {
        ac_fixed<W,I,S> r = this->operator>>( op2).val_;
        *this = r;
        return *this;
    }

    // Relational ---------------------------------------------------------------
	#pragma calypto_flag AC_FIXED_REL_EQ
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    bool operator == ( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };
        if(F == F2)
            return val_ ==  op2.val_;
        else if(F > F2)
            return (val_ == (op2.val_ <<(F-F2)));
        else
            return ((val_ << (F2-F)) == op2.val_);
    }
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    bool operator != ( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        return !(this->operator==(op2));
    }
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    bool operator < ( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };
        if(F == F2)
            return this->val_ < op2.val_;
        else if(F > F2)
            return (this->val_ < (op2.val_ << (F - F2)));
        else
            return ((this->val_ << (F2 - F)) < op2.val_);
    }
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    bool operator >= ( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const {
        return !(this->operator<(op2));
    }

    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    bool operator > ( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        enum { F=W-I, F2=W2-I2 };
        if(F == F2)
            return this->val_ > op2.val_;
        else if(F > F2)
            return (this->val_ > (op2.val_ << (F - F2)));
        else
            return ((this->val_ << (F2 - F) ) > op2.val_);
    }
    template<int W2, int I2, bool S2, ac_q_mode Q2, ac_o_mode O2>
    bool operator <= ( const ac_fixed<W2,I2,S2,Q2,O2> &op2) const
    {
        return (!(this->operator>(op2)));
    }

    /*ac_fixed& operator = ( double d) const
    {
        calypto_double d = d_val;
        double_to_ac_fixed_ctor(d);
        return *this;
    }*/

    #pragma calypto_flag AC_FIXED_ASSIGN_OP
    ac_fixed& operator = ( const calypto_double& d)
    {
        double_to_ac_fixed_ctor(d);
        return *this;

    }

	#pragma calypto_flag AC_FIXED_REL_EQ
    bool operator == ( const calypto_double& d) const
    {
        if(is_neg() != d < 0.0)
            return false;
        calypto_double di = d;
        ldexpr<-(I+(32-(W&31)))>(di);
        bool overflow, qb, r;
        ac_fixed<W,I,S> t;
        conv_from_fraction<(W+ 31)/32, W + !S, S>(di, &t.val_, &qb, &r, &overflow);
        if(qb || r || overflow)
            return false;
        return operator == (t);
    }

    bool operator < ( const calypto_double& d) const
    {
        if(is_neg() != d < 0.0)
            return is_neg();
        calypto_double di = d;
        ldexpr<-(I+(32-(W&31)))>(di);
        bool overflow, qb, r;
        ac_fixed<W,I,S> t;
        conv_from_fraction<(W +31)/32, W  + !S, S>(di, &t.val_, &qb, &r, &overflow);
        if(is_neg() && overflow)
            return false;
        return !is_neg() && overflow || (qb || r)
            && operator <= (t) || operator < (t);
    }

    bool operator >= ( const calypto_double& d) const
    {
        return !operator < ( d );
    }
    bool operator > ( const calypto_double& d) const
    {
        if(is_neg() != d < 0.0)
            return !is_neg();
        calypto_double di = d;
        ldexpr<-(I+(32-(W&31)))>(di);
        bool overflow, qb, r;
        ac_fixed<W,I,S> t;
        conv_from_fraction<(W +31)/32, W  + !S, S>(di, &t.val_, &qb, &r, &overflow);
        if(!is_neg() && overflow )
            return false;
        return is_neg() && overflow || operator > (t);
    }

    bool operator != ( const calypto_double& d) const
    {
        return !operator == ( d );
    }

    bool operator <= ( const calypto_double& d) const
    {
        return !operator > ( d );
    }

    // Bit and Slice Select --------------------------------------
    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_SLICE
    template<int WS, int WX, bool SX>
    inline ac_int<WS,S> slc(const ac_int<WX,SX> &index) const
    {
#ifdef CALYPTO_ASSERT
        calypto_assert(index >= 0,
                "Attempting to read slc with negative indeces");
#endif
        enum { ShiftWidth = AC_MIN(WX-SX,     // Strip off the sign bit
                                   32-SX) };  // Catapult considers at most 32 bits (minus sign)
        sc_biguint<ShiftWidth> shift = to_sc(index);

#ifdef SLEC_CPC
        unsigned uindex = (shift.to_uint()) & ((unsigned)~0 >> 1);
        return slc<WS>( uindex );
#else
        ac_fixed_dt<W,S>::dt temp = val_ >> shift;
        ac_int<WS,S> r(temp);
        return r;
#endif
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_SLICE
    template<int WS>
    inline ac_int<WS,S> slc(signed index) const
    {
#ifdef CALYPTO_ASSERT
        calypto_assert(index >= 0, "Attempting to read slc with negative indeces");
#endif
        unsigned uindex = index & ((unsigned)~0 >> 1);
#ifdef SLEC_CPC
        ac_fixed_dt<WS,S>::dt value(val_.range(uindex+(WS-1),uindex));
        ac_fixed_dt<W,S>::dt temp = value;
#else
        ac_fixed_dt<W,S>::dt temp =val_ >> uindex ;
#endif
        ac_int<WS,S> r ( temp) ;
        return r;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_SLICE
    template<int WS>
    inline ac_int<WS,S> slc(unsigned uindex) const
    {
#ifdef SLEC_CPC
        uindex = uindex & ((unsigned)~0 >> 1);
        ac_fixed_dt<WS,S>::dt value(val_.range(uindex+(WS-1),uindex));
        ac_fixed_dt<W,S>::dt temp = value;
#else
        ac_fixed_dt<W,S>::dt temp = val_ >> uindex;
#endif
        ac_int<WS,S> r(temp);
        return r;
    }


    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_SET_SLICE
    template<int W2, bool S2, int WX, bool SX>
    inline ac_fixed &set_slc(const ac_int<WX,SX> lsb,
        const ac_int<W2,S2> &slc)
    {
        return set_slc(lsb.to_uint(), slc);
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_SET_SLICE
    template<int W2, bool S2>
    inline ac_fixed &set_slc(signed lsb, const ac_int<W2,S2> &slc)
    {
        unsigned ulsb = lsb & ((unsigned)~0 >> 1);
        return set_slc(ulsb, slc);
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_SET_SLICE
    template<int W2, bool S2>
    inline ac_fixed &set_slc(unsigned ulsb, const ac_int<W2,S2> &slc)
    {
        // Warning: cat2slec will mismatch if OOB writes occur here.
#ifdef CALYPTO_ASSERT
        calypto_assert((ulsb >= 0) && (ulsb + W2 <= W), "Attempting an out-of-bounds write");
#endif
        this->val_.range (ulsb+(W2-1), ulsb) = to_sc(slc);
        return *this;
    }

    class ac_bitref {
        ac_fixed* d_bv;
        unsigned d_index;
    public:
        ac_bitref()  : d_bv(NULL), d_index(0) {}
        ac_bitref( ac_fixed *bv, unsigned index=0 )
            : d_bv(bv), d_index(index)
        {}
        operator bool () const
        {
            if( (d_index < W) )
            {
                return (d_index < W) ?(d_bv->val_[ d_index]!=0) : 0;
                //return (d_index < W) ?(((d_bv.val_ >> d_index ) & 1)!=0) : 0;
            }
            return false;
        }

        #pragma calypto_flag CALLREPLACEPRAGMA
        inline ac_bitref operator = ( int val )
        {
            if ((d_bv != NULL) && (d_index < W))
            {
                d_bv->val_[d_index] = val&1;
            }
            return *this;

        }
        template<int W2, bool S2>
        inline ac_bitref operator = ( const ac_int<W2,S2> &val )
        {
            return operator =(val.to_int());
        }
        inline ac_bitref operator = ( const ac_bitref &val )
        {
            return operator =((int) (bool) val);
        }
    };

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_BIT_ACCESS
    ac_bitref operator [] ( unsigned int uindex)
    {
        ac_bitref bvh( this, uindex );
        return bvh;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_BIT_ACCESS
    ac_bitref operator [] ( int index)
    {
        unsigned uindex = index & ((unsigned)~0 >> 1);
        ac_bitref bvh( this, uindex );
        return bvh;
    }
    #pragma calypto_flag AC_FIXED_BIT_ACCESS
    template<int W2, bool S2>
    ac_bitref operator [] ( const ac_int<W2,S2> &index)
    {
        ac_int<W2-S2,false> uindex = index;
        ac_bitref bvh( this, uindex.to_uint() );
        return bvh;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_BIT_ACCESS
    bool operator [] ( unsigned int uindex) const
    {
#ifdef CALYPTO_ASSERT
        calypto_assert(uindex < W, "Attempting to read bit beyond MSB");
#endif
        //return (uindex < W) ?(((val_ >> uindex ) & 1)!=0) : 0;
        return (uindex < W) ?((val_[uindex] )!=0) : 0;
    }

    #pragma calypto_flag CALLREPLACEPRAGMA
    #pragma calypto_flag AC_FIXED_BIT_ACCESS
    bool operator [] ( int index) const
    {
#ifdef CALYPTO_ASSERT
        calypto_assert(index >= 0,
            "Attempting to read bit with negative index");
        calypto_assert(index < W, "Attempting to read bit beyond MSB");
#endif
        unsigned uindex = index & ((unsigned)~0 >> 1);
        return (uindex < W) ? ((val_[uindex])!=0) : 0;
        //return (uindex < W) ? (((val_ >> uindex ) & 1)!=0) : 0;

    }
    #pragma calypto_flag AC_FIXED_BIT_ACCESS
    template<int W2, bool S2>
    bool operator [] ( const ac_int<W2,S2> &index) const
    {
#ifdef CALYPTO_ASSERT
        calypto_assert(index >= 0, "Attempting to read bit with negative index");
        calypto_assert(index < W, "Attempting to read bit beyond MSB");
#endif
        ac_int<W2-S2,false> uindex = index;

        return (uindex < W) ? ((val_[uindex])!=0) : 0;
        //return (uindex < W) ? (((val_ >> uindex ) & 1)!=0) : 0;
    }

    inline std::string to_string(ac_base_mode base_rep, bool sign_mag = false) const {
    }

    unsigned leading_sign()
    {
        unsigned ls = 0;
        bool msb = this->operator[](W-1);

        if ( !S && msb )
        {
            return ls;
        }

        #pragma unroll yes
        for ( int i = W-2; i >= 0; --i )
        {
            if ( this->operator[](i) != msb )
            {
                break;
            }

            ++ls;
        }

        return ls + !S;
    }

    unsigned leading_sign( bool& all_sign )
    {
        unsigned ls = leading_sign();
        all_sign = ( ls == W - S );
        return ls;
    }

    inline void bit_fill_hex(const char *str) 
    {
        // Zero Pads if str is too short, throws ms bits away if str is too long
        // Asserts if anything other than 0-9a-fA-F is encountered
        ac_int<W,S> x;
        x.bit_fill_hex(str);
        set_slc(0, x);
    }
    
    template<int N>
    inline void bit_fill(const int (&ivec)[N], bool bigendian=true)
    {
        // bit_fill from integer vector
        //   if W > N*32, missing most significant bits are zeroed
        //   if W < N*32, additional bits in ivec are ignored (no overflow checking)
        //
        // Example:  
        //   ac_fixed<80,40,false> x;    int vec[] = { 0xffffa987, 0x6543210f, 0xedcba987 };
        //   x.bit_fill(vec);   // vec[0] fill bits 79-64 
        ac_int<W,S> x;
        x.bit_fill(ivec, bigendian);
        set_slc(0, x);
    }

};

namespace ac {
  template<typename T>
  struct ac_fixed_represent {
    enum { t_w = ac_private::c_type_params<T>::W, t_i = t_w, t_s = ac_private::c_type_params<T>::S };
    typedef ac_fixed<t_w,t_i,t_s> type;
  };
  template<> struct ac_fixed_represent<float> {};
  template<> struct ac_fixed_represent<double> {};
  template<int W, bool S>
  struct ac_fixed_represent< ac_int<W,S> > {
    typedef ac_fixed<W,W,S> type;
  };
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
  struct ac_fixed_represent< ac_fixed<W,I,S,Q,O> > {
    typedef ac_fixed<W,I,S,Q,O> type;
  };
}

namespace ac_private {
  // with T == ac_fixed
  template<int W2, int I2, bool S2>
  struct rt_ac_fixed_T< ac_fixed<W2,I2,S2> > {
    template<int W, int I, bool S>
    struct op1 {
      typedef typename ac_fixed<W,I,S>::template rt<W2,I2,S2>::mult mult;
      typedef typename ac_fixed<W,I,S>::template rt<W2,I2,S2>::plus plus;
      typedef typename ac_fixed<W,I,S>::template rt<W2,I2,S2>::minus minus;
      typedef typename ac_fixed<W2,I2,S2>::template rt<W,I,S>::minus minus2;
      typedef typename ac_fixed<W,I,S>::template rt<W2,I2,S2>::logic logic;
      typedef typename ac_fixed<W,I,S>::template rt<W2,I2,S2>::div div;
      typedef typename ac_fixed<W2,I2,S2>::template rt<W,I,S>::div div2;
    };
  };
  // with T == ac_int
  template<int W2, bool S2>
  struct rt_ac_fixed_T< ac_int<W2,S2> > {
    template<int W, int I, bool S>
    struct op1 {
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::mult mult;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::plus plus;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::minus minus;
      typedef typename ac_fixed<W2,W2,S2>::template rt<W,I,S>::minus minus2;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::logic logic;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::div div;
      typedef typename ac_fixed<W2,W2,S2>::template rt<W,I,S>::div div2;
    };
  };

  template<typename T>
  struct rt_ac_fixed_T< c_type<T> > {
    enum { W2 = c_type_params<T>::W, I2 = W2, S2 = c_type_params<T>::S };
    template<int W, int I, bool S>
    struct op1 {
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::mult mult;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::plus plus;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::minus minus;
      typedef typename ac_fixed<W2,W2,S2>::template rt<W,I,S>::minus minus2;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::logic logic;
      typedef typename ac_fixed<W,I,S>::template rt<W2,W2,S2>::div div;
      typedef typename ac_fixed<W2,W2,S2>::template rt<W,I,S>::div div2;
    };
  };
}

template<> inline ac_fixed<1,1,true,AC_TRN,AC_WRAP>::ac_fixed( bool b ) { val_ = b ? -1 : 0; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( bool b ) { val_ = b; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( signed char b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned char b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( signed short b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned short b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( signed int b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned int b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( signed long b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned long b ) { val_ = b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( Ulong b ) { val_ = (int) b&1; }

template<> inline ac_fixed<1,1,false,AC_TRN,AC_WRAP>::ac_fixed( Slong b ) { val_ = (int) b&1; }

template<> inline ac_fixed<8,8,true,AC_TRN,AC_WRAP>::ac_fixed( bool b ) { val_ = b; }

template<> inline ac_fixed<8,8,false,AC_TRN,AC_WRAP>::ac_fixed( bool b ) { val_ = b; }

template<> inline ac_fixed<8,8,true,AC_TRN,AC_WRAP>::ac_fixed( signed char b ) { val_ = b; }

template<> inline ac_fixed<8,8,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned char b ) { val_ = b; }

template<> inline ac_fixed<8,8,true,AC_TRN,AC_WRAP>::ac_fixed( unsigned char b ) { val_ = (signed char) b; }

template<> inline ac_fixed<8,8,false,AC_TRN,AC_WRAP>::ac_fixed( signed char b ) { val_ = (unsigned char) b; }

template<> inline ac_fixed<16,16,true,AC_TRN,AC_WRAP>::ac_fixed( bool b ) { val_ = b; }

template<> inline ac_fixed<16,16,false,AC_TRN,AC_WRAP>::ac_fixed( bool b ) { val_ = b; }

template<> inline ac_fixed<16,16,true,AC_TRN,AC_WRAP>::ac_fixed( signed char b ) { val_ = b; }

template<> inline ac_fixed<16,16,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned char b ) { val_ = b; }

template<> inline ac_fixed<16,16,true,AC_TRN,AC_WRAP>::ac_fixed( unsigned char b ) { val_ = b; }

template<> inline ac_fixed<16,16,false,AC_TRN,AC_WRAP>::ac_fixed( signed char b ) { val_ = (unsigned short) b; }

template<> inline ac_fixed<16,16,true,AC_TRN,AC_WRAP>::ac_fixed( signed short b ) { val_ = b; }

template<> inline ac_fixed<16,16,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned short b ) { val_ = b; }

template<> inline ac_fixed<16,16,true,AC_TRN,AC_WRAP>::ac_fixed( unsigned short b ) { val_ = (signed short) b; }

template<> inline ac_fixed<16,16,false,AC_TRN,AC_WRAP>::ac_fixed( signed short b ) { val_ = (unsigned short) b; }

template<> inline ac_fixed<32,32,true,AC_TRN,AC_WRAP>::ac_fixed( signed int b ) { val_ = b; }

template<> inline ac_fixed<32,32,true,AC_TRN,AC_WRAP>::ac_fixed( unsigned int b ) { val_ = b; }

template<> inline ac_fixed<32,32,false,AC_TRN,AC_WRAP>::ac_fixed( signed int b ) { val_ = (unsigned int) b; }

template<> inline ac_fixed<32,32,false,AC_TRN,AC_WRAP>::ac_fixed( unsigned int b ) { val_ = b; }

template<> inline ac_fixed<32,32,true,AC_TRN,AC_WRAP>::ac_fixed( Slong b ) { val_ = (int) b; }

template<> inline ac_fixed<32,32,true,AC_TRN,AC_WRAP>::ac_fixed( Ulong b ) { val_ = (int) b; }

template<> inline ac_fixed<32,32,false,AC_TRN,AC_WRAP>::ac_fixed( Slong b ) { val_ = (unsigned int) b; }

template<> inline ac_fixed<32,32,false,AC_TRN,AC_WRAP>::ac_fixed( Ulong b ) { val_ = (unsigned int) b; }

template<> inline ac_fixed<64,64,true,AC_TRN,AC_WRAP>::ac_fixed( Slong b ) { val_ =  b;  }

template<> inline ac_fixed<64,64,true,AC_TRN,AC_WRAP>::ac_fixed( Ulong b ) { val_ =  b; }

template<> inline ac_fixed<64,64,false,AC_TRN,AC_WRAP>::ac_fixed( Slong b ) { val_ =  (unsigned long long) b; }

template<> inline ac_fixed<64,64,false,AC_TRN,AC_WRAP>::ac_fixed( Ulong b ) { val_ =  (unsigned long long) b; }

// Stream --------------------------------------------------------------------

template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline std::ostream& operator << (std::ostream &os, const ac_fixed<W,I,S,Q,O> &x) {
  return os;
}

 // Binary Operators with Integers --------------------------------------------

#define FX_BIN_OP_WITH_INT_2I(BIN_OP, C_TYPE, WI, SI, RTYPE)  \
  _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline typename ac_fixed<W,I,S>::template rt<WI,WI,SI>::RTYPE operator BIN_OP ( const ac_fixed<W,I,S,Q,O> &op, C_TYPE i_op) {  \
    return op.operator BIN_OP (ac_int<WI,SI>(i_op));  \
  }

#define FX_BIN_OP_WITH_INT(BIN_OP, C_TYPE, WI, SI, RTYPE)  \
    FX_BIN_OP_WITH_INT1(BIN_OP, C_TYPE, WI, SI, RTYPE)  \
    FX_BIN_OP_WITH_INT2(BIN_OP, C_TYPE, WI, SI, RTYPE)


#define FX_BIN_OP_WITH_INT_DIV(BIN_OP, C_TYPE, WI, SI, RTYPE)  \
   _Pragma("calypto_flag AC_FIXED_DIVIDE") \
   FX_BIN_OP_WITH_INT1(BIN_OP, C_TYPE, WI, SI, RTYPE)  \
   _Pragma("calypto_flag AC_FIXED_DIVIDE") \
   FX_BIN_OP_WITH_INT2(BIN_OP, C_TYPE, WI, SI, RTYPE)

#define FX_BIN_OP_WITH_INT1(BIN_OP, C_TYPE, WI, SI, RTYPE)  \
  _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline typename ac_fixed<WI,WI,SI>::template rt<W,I,S>::RTYPE operator BIN_OP ( C_TYPE i_op, const ac_fixed<W,I,S,Q,O> &op) {  \
    return ac_fixed<WI,WI,SI>(i_op).operator BIN_OP (op);  \
  } 

#define  FX_BIN_OP_WITH_INT2(BIN_OP, C_TYPE, WI, SI, RTYPE)  \
  _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline typename ac_fixed<W,I,S>::template rt<WI,WI,SI>::RTYPE operator BIN_OP ( const ac_fixed<W,I,S,Q,O> &op, C_TYPE i_op) {  \
    return op.operator BIN_OP (ac_fixed<WI,WI,SI>(i_op));  \
  }

#define FX_REL_OP_WITH_INT(REL_OP, C_TYPE, W2, S2)  \
  _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline bool operator REL_OP ( const ac_fixed<W,I,S,Q,O> &op, C_TYPE op2) {  \
    return op.operator REL_OP (ac_fixed<W2,W2,S2>(op2));  \
  }  \
  _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline bool operator REL_OP ( C_TYPE op2, const ac_fixed<W,I,S,Q,O> &op) {  \
    return ac_fixed<W2,W2,S2>(op2).operator REL_OP (op);  \
  }

#define FX_ASSIGN_OP_WITH_INT_2(ASSIGN_OP, C_TYPE, W2, S2)  \
  _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline ac_fixed<W,I,S,Q,O> &operator ASSIGN_OP ( ac_fixed<W,I,S,Q,O> &op, C_TYPE op2) {  \
    return op.operator ASSIGN_OP (ac_fixed<W2,W2,S2>(op2));  \
  }

#define FX_ASSIGN_OP_WITH_INT_2I(ASSIGN_OP, C_TYPE, W2, S2)  \
  _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline ac_fixed<W,I,S,Q,O> &operator ASSIGN_OP ( ac_fixed<W,I,S,Q,O> &op, C_TYPE op2) {  \
    return op.operator ASSIGN_OP (ac_int<W2,S2>(op2));  \
  }

#define FX_OPS_WITH_INT(C_TYPE, WI, SI) \
  FX_BIN_OP_WITH_INT(*, C_TYPE, WI, SI, mult) \
  FX_BIN_OP_WITH_INT(+, C_TYPE, WI, SI, plus) \
  FX_BIN_OP_WITH_INT(-, C_TYPE, WI, SI, minus) \
  FX_BIN_OP_WITH_INT_DIV(/, C_TYPE, WI, SI, div) \
  FX_BIN_OP_WITH_INT_2I(>>, C_TYPE, WI, SI, arg1) \
  FX_BIN_OP_WITH_INT_2I(<<, C_TYPE, WI, SI, arg1) \
  FX_BIN_OP_WITH_INT(&, C_TYPE, WI, SI, logic) \
  FX_BIN_OP_WITH_INT(|, C_TYPE, WI, SI, logic) \
  FX_BIN_OP_WITH_INT(^, C_TYPE, WI, SI, logic) \
  \
  FX_REL_OP_WITH_INT(==, C_TYPE, WI, SI) \
  FX_REL_OP_WITH_INT(!=, C_TYPE, WI, SI) \
  FX_REL_OP_WITH_INT(>, C_TYPE, WI, SI) \
  FX_REL_OP_WITH_INT(>=, C_TYPE, WI, SI) \
  FX_REL_OP_WITH_INT(<, C_TYPE, WI, SI) \
  FX_REL_OP_WITH_INT(<=, C_TYPE, WI, SI) \
  \
  FX_ASSIGN_OP_WITH_INT_2(+=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2(-=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2(*=, C_TYPE, WI, SI) \
   _Pragma("calypto_flag AC_FIXED_DIVEQUALS") \
  FX_ASSIGN_OP_WITH_INT_2(/=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2(%=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2I(>>=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2I(<<=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2(&=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2(|=, C_TYPE, WI, SI) \
  FX_ASSIGN_OP_WITH_INT_2(^=, C_TYPE, WI, SI)


 // Operators with Double --------------------------------------------

#define FX_REL_OP_WITH_DOUBLE(REL_OP, C_TYPE, W2,I2, S2)  \
    _Pragma("calypto_flag CALYPTO_FLOAT_CONST_OP") \
    _Pragma("calypto_flag CALLREPLACEPRAGMA")  \
    _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline bool operator REL_OP ( const ac_fixed<W,I,S,Q,O> &op, C_TYPE op2) {  \
    return op.operator REL_OP (ac_fixed<W2,I2,S2>(op2));  \
  }  \
    _Pragma("calypto_flag CALYPTO_FLOAT_CONST_OP") \
    _Pragma("calypto_flag CALLREPLACEPRAGMA")  \
    _Pragma("calypto_flag AC_FIXED_WITH_CTYPES_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O> \
  inline bool operator REL_OP ( C_TYPE op2, const ac_fixed<W,I,S,Q,O> &op) {  \
    return ac_fixed<W2,I2,S2>(op2).operator REL_OP (op);  \
  }



  #define FX_OPS_WITH_DOUBLE(C_TYPE, WI ,I1, SI) \
  FX_REL_OP_WITH_DOUBLE(==, C_TYPE, WI,I1, SI) \
  FX_REL_OP_WITH_DOUBLE(!=, C_TYPE, WI,I1, SI) \
  FX_REL_OP_WITH_DOUBLE(>, C_TYPE, WI, I1,SI) \
  FX_REL_OP_WITH_DOUBLE(>=, C_TYPE, WI,I1, SI) \
  FX_REL_OP_WITH_DOUBLE(<, C_TYPE, WI, I1,SI) \
  FX_REL_OP_WITH_DOUBLE(<=, C_TYPE, WI, I1,SI) \



// --------------------------------------- End of Macros for Binary Operators with C Integers

namespace ac {
  namespace ops_with_other_types {
    // Binary Operators with C Integers --------------------------------------------
    FX_OPS_WITH_INT(bool, 1, false)
    FX_OPS_WITH_INT(char, 8, true)
    FX_OPS_WITH_INT(signed char, 8, true)
    FX_OPS_WITH_INT(unsigned char, 8, false)
    FX_OPS_WITH_INT(short, 16, true)
    FX_OPS_WITH_INT(unsigned short, 16, false)
    FX_OPS_WITH_INT(int, 32, true)
    FX_OPS_WITH_INT(unsigned int, 32, false)
    FX_OPS_WITH_INT(long, 32, true)
    FX_OPS_WITH_INT(unsigned long, 32, false)
    FX_OPS_WITH_INT(Slong, 64, true)
    FX_OPS_WITH_INT(Ulong, 64, false)
    FX_OPS_WITH_DOUBLE(double, 128, 64, true)
    // -------------------------------------- End of Binary Operators with Integers
  }  // ops_with_other_types namespace
} // ac namespace


// Macros for Binary Operators with ac_int --------------------------------------------

#define FX_BIN_OP_WITH_AC_INT_1(BIN_OP, RTYPE)  \
  _Pragma("calypto_flag AC_FIXED_WITH_AC_INT_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int WI, bool SI> \
  inline typename ac_fixed<WI,WI,SI>::template rt<W,I,S>::RTYPE operator BIN_OP ( const ac_int<WI,SI> &i_op, const ac_fixed<W,I,S,Q,O> &op) {  \
    return ac_fixed<WI,WI,SI>(i_op).operator BIN_OP (op);  \
  }

#define FX_BIN_OP_WITH_AC_INT_2(BIN_OP, RTYPE)  \
  _Pragma("calypto_flag AC_FIXED_WITH_AC_INT_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int WI, bool SI> \
  inline typename ac_fixed<W,I,S>::template rt<WI,WI,SI>::RTYPE operator BIN_OP ( const ac_fixed<W,I,S,Q,O> &op, const ac_int<WI,SI> &i_op) {  \
    return op.operator BIN_OP (ac_fixed<WI,WI,SI>(i_op));  \
  }

#define FX_BIN_OP_WITH_AC_INT(BIN_OP, RTYPE)  \
  FX_BIN_OP_WITH_AC_INT_1(BIN_OP, RTYPE) \
  FX_BIN_OP_WITH_AC_INT_2(BIN_OP, RTYPE)

#define FX_BIN_OP_WITH_AC_INT_DIV(BIN_OP, RTYPE)  \
  _Pragma("calypto_flag AC_FIXED_DIVIDE") \
  FX_BIN_OP_WITH_AC_INT_1(BIN_OP, RTYPE) \
  _Pragma("calypto_flag AC_FIXED_DIVIDE") \
  FX_BIN_OP_WITH_AC_INT_2(BIN_OP, RTYPE)

#define FX_REL_OP_WITH_AC_INT(REL_OP)  \
  _Pragma("calypto_flag AC_FIXED_WITH_AC_INT_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int WI, bool SI> \
  inline bool operator REL_OP ( const ac_fixed<W,I,S,Q,O> &op, const ac_int<WI,SI> &op2) {  \
    return op.operator REL_OP (ac_fixed<WI,WI,SI>(op2));  \
  }  \
  _Pragma("calypto_flag AC_FIXED_WITH_AC_INT_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int WI, bool SI> \
  inline bool operator REL_OP ( ac_int<WI,SI> &op2, const ac_fixed<W,I,S,Q,O> &op) {  \
    return ac_fixed<WI,WI,SI>(op2).operator REL_OP (op);  \
  }

#define FX_ASSIGN_OP_WITH_AC_INT(ASSIGN_OP) \
    FX_ASSIGN_OP_WITH_AC_INT1(ASSIGN_OP) \
    FX_ASSIGN_OP_WITH_AC_INT2(ASSIGN_OP)

#define FX_ASSIGN_OP_WITH_AC_INT_DIVEQ(ASSIGN_OP) \
    _Pragma("calypto_flag AC_FIXED_DIVEQUALS") \
    FX_ASSIGN_OP_WITH_AC_INT1(ASSIGN_OP) \
    _Pragma("calypto_flag AC_FIXED_DIVEQUALS") \
    FX_ASSIGN_OP_WITH_AC_INT2(ASSIGN_OP) 

#define FX_ASSIGN_OP_WITH_AC_INT1(ASSIGN_OP)  \
  _Pragma("calypto_flag AC_FIXED_WITH_AC_INT_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int WI, bool SI> \
  inline ac_fixed<W,I,S,Q,O> &operator ASSIGN_OP ( ac_fixed<W,I,S,Q,O> &op, const ac_int<WI,SI> &op2) {  \
    return op.operator ASSIGN_OP (ac_fixed<WI,WI,SI>(op2));  \
  }  

#define FX_ASSIGN_OP_WITH_AC_INT2(ASSIGN_OP)  \
  _Pragma("calypto_flag AC_FIXED_WITH_AC_INT_BINOP") \
  template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O, int WI, bool SI> \
  inline ac_int<WI,SI> &operator ASSIGN_OP ( ac_int<WI,SI> &op, const ac_fixed<W,I,S,Q,O> &op2) {  \
    return op.operator ASSIGN_OP (op2.to_ac_int());  \
  }

namespace ac {
  namespace ops_with_other_types {
    // Binary Operators with ac_int --------------------------------------------

FX_BIN_OP_WITH_AC_INT(*, mult)
FX_BIN_OP_WITH_AC_INT(+, plus)
FX_BIN_OP_WITH_AC_INT(-, minus)
FX_BIN_OP_WITH_AC_INT_DIV(/, div)
FX_BIN_OP_WITH_AC_INT(&, logic)
FX_BIN_OP_WITH_AC_INT(|, logic)
FX_BIN_OP_WITH_AC_INT(^, logic)

FX_REL_OP_WITH_AC_INT(==)
FX_REL_OP_WITH_AC_INT(!=)
FX_REL_OP_WITH_AC_INT(>)
FX_REL_OP_WITH_AC_INT(>=)
FX_REL_OP_WITH_AC_INT(<)
FX_REL_OP_WITH_AC_INT(<=)

FX_ASSIGN_OP_WITH_AC_INT(+=)
FX_ASSIGN_OP_WITH_AC_INT(-=)
FX_ASSIGN_OP_WITH_AC_INT(*=)
FX_ASSIGN_OP_WITH_AC_INT_DIVEQ(/=)
FX_ASSIGN_OP_WITH_AC_INT(%=)
FX_ASSIGN_OP_WITH_AC_INT(&=)
FX_ASSIGN_OP_WITH_AC_INT(|=)
FX_ASSIGN_OP_WITH_AC_INT(^=)

// Relational Operators with double
#pragma calypto_flag AC_FIXED_REL_EQ
#pragma calypto_flag AC_FIXED_WITH_CTYPES_BINOP
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline bool operator == ( calypto_double const& op, const ac_fixed<W,I,S,Q,O> &op2) {
  return op2.operator == (op);
}
#pragma calypto_flag AC_FIXED_WITH_CTYPES_BINOP
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline bool operator != ( calypto_double const& op, const ac_fixed<W,I,S,Q,O> &op2) {
  return op2.operator != (op);
}
#pragma calypto_flag AC_FIXED_WITH_CTYPES_BINOP
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline bool operator > ( calypto_double const& op, const ac_fixed<W,I,S,Q,O> &op2) {
  return op2.operator < (op);
}
#pragma calypto_flag AC_FIXED_WITH_CTYPES_BINOP
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline bool operator < ( calypto_double const& op, const ac_fixed<W,I,S,Q,O> &op2) {
  return op2.operator > (op);
}
#pragma calypto_flag AC_FIXED_WITH_CTYPES_BINOP
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline bool operator <= ( calypto_double const& op, const ac_fixed<W,I,S,Q,O> &op2) {
  return op2.operator >= (op);
}
#pragma calypto_flag AC_FIXED_WITH_CTYPES_BINOP
template<int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline bool operator >= ( calypto_double const& op, const ac_fixed<W,I,S,Q,O> &op2) {
  return op2.operator <= (op);
}

}  // ops_with_other_types namespace
} // ac namespace

using namespace ac::ops_with_other_types;
template<ac_special_val V, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
inline ac_fixed<W,I,S,Q,O> value(ac_fixed<W,I,S,Q,O>) {
    ac_fixed<W,I,S> r;
    r.set_val<V>();
    return r;
}

namespace ac {
// PUBLIC FUNCTIONS
// function to initialize (or uninitialize) arrays
#pragma calypto_flag AC_INIT_ARRAY
template<ac_special_val V, int W, int I, bool S, ac_q_mode Q, ac_o_mode O>
  inline bool init_array(ac_fixed<W,I,S,Q,O> *a, int n) {
    ac_fixed<W,I,S> t = value<V>(*a);
    for(int i=0; i < n; i++)
      a[i] = t;
    return true;
  }

  #pragma calypto_flag API_NOT_SUPPORTED
  inline ac_fixed<54,1,true,AC_TRN, AC_WRAP> frexp_d( double d, ac_int<11,true>& exp )
  {
      ac_fixed<54,1,true,AC_TRN, AC_WRAP> ret;
      return ret;
  }

  #pragma calypto_flag API_NOT_SUPPORTED
  inline ac_fixed<25,1,true,AC_TRN, AC_WRAP> frexp_f( float f, ac_int<8,true> &exp )
  {
      ac_fixed<25,1,true,AC_TRN, AC_WRAP> ret;
      return ret;
  }

  #pragma calypto_flag API_NOT_SUPPORTED
  inline ac_fixed<53,1,false,AC_TRN, AC_WRAP> frexp_sm_d(double d, ac_int<11,true> &exp, bool &sign)
  {
      ac_fixed<53,1,false,AC_TRN, AC_WRAP> ret;
      return ret;
  }

  #pragma calypto_flag API_NOT_SUPPORTED
  inline ac_fixed<24,1,false,AC_TRN, AC_WRAP> frexp_sm_f(float f, ac_int<8,true> &exp, bool &sign)
  {
      ac_fixed<24,1,false,AC_TRN, AC_WRAP>ret;
      return ret;
  }

}//namespace ac
#endif
