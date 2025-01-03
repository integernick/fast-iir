#pragma once

#include "BiquadBlock.h"

#include <array>

namespace tiny_iir {
    template<unsigned int ORDER,
            typename T = double,
            class BiquadBlockDirectForm = BiquadBlockDF1<T>>
    class CascadeFilter {
    public:
        static constexpr unsigned int NUMBER_OF_BIQUAD_BLOCKS = (ORDER + 1) / 2;
        static constexpr unsigned int COEFFICIENTS_PER_BIQUAD_BLOCK = 5;
        static constexpr unsigned int NUMBER_OF_COEFFICIENTS
            = NUMBER_OF_BIQUAD_BLOCKS * COEFFICIENTS_PER_BIQUAD_BLOCK;

        explicit CascadeFilter(const T coefficients[NUMBER_OF_COEFFICIENTS]) {
            load_coefficients(coefficients);
            update_coefficients();
        }

        [[nodiscard]] T *get_coefficients() const {
            return _coefficients;
        }

        void load_coefficients(std::array<T, NUMBER_OF_COEFFICIENTS> coefficients) {
            memcpy(_coefficients, coefficients, sizeof(T) * NUMBER_OF_COEFFICIENTS);
        }

        void update_coefficients() {
            for (unsigned int i = 0; i < NUMBER_OF_BIQUAD_BLOCKS; ++i) {
                _biquad_blocks[i].set_coefficients(_coefficients + i * COEFFICIENTS_PER_BIQUAD_BLOCK);
            }
        }

        [[nodiscard]] T get_gain() const {
            return _gain;
        }

        void set_gain(T gain) {
            _gain = gain;
        }

        T process(T x) {
            x *= _gain;
            for (auto &biquad_block : _biquad_blocks) {
                x = biquad_block.process(x);
            }
            x *= _gain;
            return x;
        }

    private:
        T _gain;
        T _coefficients[NUMBER_OF_COEFFICIENTS];
        BiquadBlockDirectForm _biquad_blocks[NUMBER_OF_BIQUAD_BLOCKS];
    };
}
