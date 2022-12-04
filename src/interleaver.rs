
/// Precompute some intermediate values and lookup tables that are needed
/// when interleaving a vector to form the bytes that will be built into a `BigUint` Hilbert Index.
/// This structure maps each bit from each u32 dimension in the input to a corresponding bit and byte in the output.
/// 
/// Note: The goal of this class was to speed up the Hilbert transformation. Benchmarks show no detectable improvement.
pub struct Interleaver {
    /// Number of bytes needed by the output.
    bytes_needed : usize,

    /// Number of bits needed by the output. 
    num_bits : usize,

    /// From which input dimension should a given output bit be taken? 
    /// As i_bit increases, this goes from first dimension to last and repeats in a modular fashion.
    i_from_uint_vector : Vec<usize>,

    /// From which input bit should a given output bit be taken? 
    /// Progress from high-bit to low bit of the input.
    i_from_uint_bit : Vec<usize>,

    /// Into which output byte should a given bit be written? 
    /// In the output, we need to skip past the pad bits.
    i_to_byte_vector : Vec<usize>,

    /// Into which output bit position should a given bit be written?
    i_to_byte_bit : Vec<usize>
}

impl Interleaver {

    pub fn new(dimensions : usize, bit_depth : usize) -> Self {
        // The individual bytes in the value array will be created in big-endian order, from highest-order byte to lowest-order byte
        // as prep for creating a BigUint in untranspose.

        // Due to rounding we may need a few bits extra of padding if bit_depth * dimensions is not an even multiple of eight bits, so add 7 bits before
        // dividing by 8 (the right shift by 3). 
        let bytes_needed = (bit_depth * dimensions + 7) >> 3;
        // BigUint (unlike C# BigInteger) does not need an extra, zero byte at the end for a sign bit, as it is unsigned.

        let num_bits = dimensions * bit_depth;

        // Number of pad bits (zeroes) that will prefix the output to make it an even number of bytes.
        let pad_bits = bytes_needed * 8 - num_bits;
        let mut i_from_uint_vector : Vec<usize> = vec![0;num_bits];
        let mut i_from_uint_bit : Vec<usize> = vec![0;num_bits];
        let mut i_to_byte_vector : Vec<usize> = vec![0;num_bits];
        let mut i_to_byte_bit : Vec<usize> = vec![0;num_bits];

        for i_bit in 0 .. num_bits
        {
            i_from_uint_vector[i_bit] = i_bit % dimensions; 
            i_from_uint_bit[i_bit] = bit_depth - (i_bit / dimensions) - 1; 
            i_to_byte_vector[i_bit] = (i_bit + pad_bits) >> 3; 
            i_to_byte_bit[i_bit] = 0x7 - ((i_bit + pad_bits) & 0x7);
        }
        Interleaver {
            bytes_needed,
            num_bits,
            i_from_uint_vector,
            i_from_uint_bit,
            i_to_byte_vector,
            i_to_byte_bit
        }
    }

    /// Interleave the bits from the vector into and return an array of bytes in big-endian order.
    /// These bytes are suitable to be combined into a single BigUint.
    pub fn interleave(&self, vector : &[u32]) -> Vec<u8>
    {
        let mut byte_vector = vec![0_u8; self.bytes_needed]; 
        for i_bit in 0 .. self.num_bits
        {
            let bit : u8 = (((vector[self.i_from_uint_vector[i_bit]] >> self.i_from_uint_bit[i_bit]) & 1_u32) << self.i_to_byte_bit[i_bit]) as u8;
            byte_vector[self.i_to_byte_vector[i_bit]] |= bit;
        }
        byte_vector
    }

}

#[cfg(test)]
mod tests {
    #[allow(unused_imports)]
    use crate::transform::fast_hilbert;
    use crate::interleaver::Interleaver;

    #[test]
    fn interleave() {
        // 3-D Axes in Binary, 5-bits per number: decimal [17, 24, 6] = binary [10001, 11000, 00110]
        // Distributed to two bytes, putting one bit padding in the high-bit position:
        //    [01100100,01001100] = [100,76] in decimal
        let dimensions = 3;
        let bit_depth = 5;
        let axes : Vec<u32> = vec![17, 24, 6];
        let interleaver = Interleaver::new(dimensions, bit_depth);
        let actual = fast_hilbert::interleave_be(&axes, 5, Some(&interleaver));
        let expected : Vec<u8> = vec![100,76];
        assert_eq!(&actual, &expected, "Correct interleave result using Interleaver");
    }

}
