//! The `transform` module holds the functions which perform the forward and inverse **Hilbert Transform**.
//! A more convenient way to prepare for and perform the forward transform is to create `Points` and
//! call `Point::hilbert_transform`. 

/// Conversions from N-dimensional points (called axes, following Skilling's usage) to Hilbert curve indices and back again.
/// 
///   - Converting from an N-dimensional point to a 1-Dimensional Hilbert Index will be called the **Hilbert Transform**.
///   - Converting from the Hilbert Index back to the N-Dimensional point will be called the **Inverse Hilbert Transform**.
/// 
/// The intermediate form of the Hilbert index is an array of transposed bits. 
/// 
/// Example: 5 bits for each of n=3 coordinates.
/// 
/// 15-bit Hilbert integer (where A is the high-order bit and O is the low-order bit) called the transposed form: 
/// 
///   A B C D E F G H I J K L M N O 
/// 
/// ```ignore
/// 
///  //                                           ↑
///  //   X[0] = A D G J M                    X[2]|  7
///  //   X[1] = B E H K N        «–––––––»       | /X[1]
///  //   X[2] = C F I L O                   axes |/
///  //          high low                         0––––––––→ X[0]
/// 
/// ```
/// NOTE: This algorithm is derived from work done by John Skilling and published in "Programming the Hilbert curve".
/// (c) 2004 American Institute of Physics.
/// 
/// NOTE: These are the most important methods:
/// 
///    1. pub fn hilbert_axes(hilbert_index : &BigUint, bits : usize, dimensions : usize) -> Vec<u32>
///       Converts Hilbert Index to N-Space.
/// 
///    2. pub fn hilbert_index(hilbert_axes : &[u32], bits  usize, interleaver_option : Option<&Interleaver>) -> BigUint
///       Converts N-Space to Hilbert Index.
/// 
/// NOTE: There are four transformations in all: 
///    
///    1. N-dimensional coordinates(&[u32]) to transposed form (Vec<u8>)
///    2. Transposed form to Hilbert index (BigUint)
///    3. Hilbert index to Transposed form 
///    4. Transposed form to N-dimensional coordinates
/// 
/// NOTE: The confusing lamination, delamination, transposed form and the like are abstracted out into a more convenient 
///       API via the `point` and `point_list` modules. Use the `Point` struct for simplicity.
pub mod fast_hilbert {
    use num::BigUint;
    use crate::interleaver::Interleaver;

    /// Convert a Hilbert distance (index) stored in a BigUint into a transposed matrix, 
    /// where the bits are distributed among many integers in an array.
    /// 
    ///   - hilbert_index - Hilbert distance as a BigUint
    ///   - bits - Number of bits per point in the N-space
    ///   - dimensions - Number of dimensions of each point in the N-space
    ///   - Returns - The Hilbert distance expressed as a transposed array.
    pub fn transpose(hilbert_index : &BigUint, bits : usize, dimensions : usize) -> Vec<u32>
    {
        // First coordinate gets its high bit from the highest bit of the hilbert_index.
        // Second coordinate gets its high bit from the second bit of the index. Etc.
        uninterleave(hilbert_index, bits, dimensions, true)
    }

    /// Convert a `BigUint` into an array of bits (from low-order bit to high-order bit) 
    /// within a framework that expects a certain number of bits.
    /// 
    ///   - n - `BigUint` to convert.
    ///   - num_bits - Mandated size of Vec to be returned. Pad the high order bits with zeroes if the BigUint 
    /// required fewer bits for its representation.
    ///   - returns - Array of ones and zeroes. The first element is the low bit of the BigInteger.
    /// 
    /// Unlike the C# version, there is no sign bit at the end.
    pub fn unpack_big_integer(n : &BigUint, num_bits : usize) -> Vec<u8>
    {
        let mut unpacked = n.to_radix_le(2);
        let unpacked_len = unpacked.len();
        // Pad the end with zeroes to ensure that we have at least num_bits.
        for _ in unpacked_len .. num_bits {
            unpacked.push(0)
        }
        unpacked
    }

    /// Spread the bits of the BigUint among coordinates for the given number of dimensions in a round-robin fashion,
    /// delaminating the big integer. 
    /// 
    /// If reverse_order is true (the normal situation) then the first (lowest order) bit of the `BigUint`
    /// is assigned to the low-order bit position of the last dimension,
    /// the next bit is assigned to the low-bit of the previous dimension, etc, until the last bit is assigned to the high-order
    /// bit position of the first dimension. This is because unpack_big_integer returns the bits with the low-order bit first. 
    /// 
    /// **Example**:
    /// 
    /// For 3 dimensions and 5 bits deep, the decimal number 25676 can be uninterleaved to the vector [17,24,6]. 
    /// Here is how: 
    /// 
    ///   -  0110010001001100 binary = 25676 in decimal
    ///   -  .1..0..0..0..1.. => 10001 = 17
    ///   -  ..1..1..0..0..0. => 11000 = 24
    ///   -  ...0..0..1..1..0 => 00110 = 6
    /// 
    /// The high-order bit is a pad bit, because we only need 15 bits, and two bytes hold 16.
    pub fn uninterleave(gray_code : &BigUint, bit_depth : usize, dimensions : usize, reverse_order : bool) -> Vec<u32>
    {
        //TODO: There must be a more efficient way to delaminate, but I can't figure it out.
        let mut axes_vector = vec![0_u32; dimensions];
        let num_bits = dimensions * bit_depth;
        let bits = unpack_big_integer(gray_code, num_bits);
        let mut bit_index = 0;
        let start_dimension : i32 = if reverse_order { dimensions as i32 - 1 } else { 0 };
        let stop_dimension : i32 = if reverse_order { -1 } else { dimensions as i32 };
        let dim_increment : i32 = if reverse_order { -1 } else { 1 };
        for bit_number in 0 .. bit_depth
        {
            let mut dimension = start_dimension;
            while dimension != stop_dimension {
                axes_vector[dimension as usize] = axes_vector[dimension as usize] | ((bits[bit_index] as u32) << bit_number);
                bit_index += 1;
                dimension += dim_increment;
            }
        }
        axes_vector
    }

    /// Convert a transposed Hilbert index back into a Big Integer (`BigUint`) index.
    /// 
    /// Assume that the number of dimensions of the point in N-space is the length of the `transposed_index` array.
    /// 
    ///   - `transposed_index` - A Hilbert index in transposed form.
    ///   - `bits` - Number of bits used to represent each coordinate.
    ///   - `interleaver_option` - Optional Interleaver to optimize the transform.
    ///   - `returns` - The Hilbert index (or distance) expressed as a `BigUint`.
    pub fn untranspose(transposed_index : &[u32], bits : usize, interleaver_option : Option<&Interleaver>) -> BigUint
    {
        // The high bit of the first coordinate becomes the high bit of the index.
        // The high bit of the second coordinate becomes the second bit of the index.

        let interleaved_bytes : Vec<u8> = interleave_be(transposed_index, bits, interleaver_option);
        // interleave_be returns a big-endian ordered array, so build the BigUint that way
        BigUint::from_bytes_be(&interleaved_bytes)
    }

    /// Interleave the bits of an unsigned vector and generate a byte array in big-endian order, 
    /// acceptable for the `BigUint` constructor (`BigUint::from_bytes_be`).
    /// 
    ///  - The high-order bit from the first number in vector becomes the high-order bit of first byte in the generated byte array (after skipping over padding bits).
    ///  - The high-order bit of the next number in vector becomes the second highest-ordered bit in the first byte in the generated byte array.
    ///  - The low-order bit of the last number becomes the low order bit of the last byte in the new array.
    /// 
    ///  - `vector` - holds the (transformed) coordinates 
    ///  - `bit_depth` - is the number of bits of precision per coordinate to use
    ///  - `interleaver_option` - Optional Interleaver to optimize the transform.
    ///  - `return` bytes in an order suitable for constructing a `BigUint`. 
    pub fn interleave_be(vector : &[u32], bit_depth : usize, interleaver_option : Option<&Interleaver>) -> Vec<u8>
    {
        if let Some(interleaver) = interleaver_option {
            return interleaver.interleave(vector);
        }

        // The individual bytes in the value array will be created in big-endian order, from highest-order byte to lowest-order byte
        // as prep for creating a BigUint in untranspose.
        let dimensions = vector.len(); // Pull member access out of loop!
        // Due to rounding we may need a few bits extra of padding if bit_depth * dimensions is not an even multiple of eight bits, so add 7 bits before
        // dividing by 8 (the right shift by 3). 
        let bytes_needed = (bit_depth * dimensions + 7) >> 3;
        // BigUint (unlike C# BigInteger) does not need an extra, zero byte at the end for a sign bit, as it is unsigned.
        let mut byte_vector = vec![0_u8; bytes_needed]; 
        let num_bits = dimensions * bit_depth;
        let pad_bits = bytes_needed * 8 - num_bits;

        for i_bit in 0 .. num_bits
        {
            let i_from_uint_vector = i_bit % dimensions; // as i_bit increases, go from first dimension to last
            let i_from_uint_bit = bit_depth - (i_bit / dimensions) - 1; // Progress from high-bit to low bit of the input.
            let i_to_byte_vector = (i_bit + pad_bits) >> 3; // In the output, we need to skip past the pad bits.
            let i_to_byte_bit = 0x7 - ((i_bit + pad_bits) & 0x7); // Progress from high-bit to low-bit of the output.

            let bit : u8 = (((vector[i_from_uint_vector] >> i_from_uint_bit) & 1_u32) << i_to_byte_bit) as u8;
            byte_vector[i_to_byte_vector] |= bit;
        }
        byte_vector
    }
/*
    fn dump_x(x : &Vec<u32>, msg : &str)
    {
        println!("X = {:?} {}", x, msg);
    }
*/
    /// Convert the Hilbert index (in its transposed form) into an N-dimensional point expressed as a vector 
    /// of unsigned coordinate values.
    /// 
    /// This performs the **Inverse Hilbert transform**. 
    /// 
    ///   - transposedIndex - The Hilbert index stored in transposed form.
    ///   - bits - Number of bits per coordinate.
    ///   - returns - Coordinate vector.
    #[allow(non_snake_case)]
    pub fn hilbert_inverse_transform(transposed_index : &[u32], bits : usize) -> Vec<u32>
    {
        let mut X : Vec<u32> = transposed_index.to_vec();
        let n = X.len(); // n: Number of dimensions
        let N : u32 = 2_u32 << (bits - 1);
        let (mut P, mut Q, mut t) : (u32, u32, u32);

        // Gray decode by H ^ (H/2)
        t = X[n - 1] >> 1;
        // Corrected error in paper which had i >= 0 for loop termination instead of i > 0 leading to negative array index.
        for i in (1..n).rev() {
            X[i] ^= X[i - 1];
        }
        X[0] ^= t;

        // Undo excess work
        Q = 2;
        while Q != N {
            P = Q - 1;
            for i in (0..n).rev() {
                if (X[i] & Q) != 0_u32 {
                    X[0] ^= P; // invert
                }
                else
                {
                    t = (X[0] ^ X[i]) & P;
                    X[0] ^= t;
                    X[i] ^= t;
                }
            }
            Q <<= 1;
        } // exchange
        X
    }

    /// Convert a Hilbert index (the distance from the origin along the Hilbert curve represented as a BigUint)
    /// into normal N-space coordinates.
    /// 
    /// Note: This is the most caller-friendly transformation from Hilbert to N-space and performs two transformations,
    /// from the BigUint to transposed form, then from transposed form to coordinates (via the **inverse Hilbert transform**). 
    /// 
    ///   - hilbertIndex - 1-dimensional Distance from origin along the Hilbert curve.
    ///   - bits - Number of bits used to encode each dimension.
    ///   - dimensions - Number of dimensions in N-space.
    ///   - returns - N-dimensional coordinate.
    pub fn hilbert_axes(hilbert_index : &BigUint, bits : usize, dimensions : usize) -> Vec<u32>
    {
        let transposed = transpose(hilbert_index, bits, dimensions);
        hilbert_inverse_transform(&transposed, bits)
    }


    /// <summary>
    /// Given the axes (coordinates) of a point in N-Dimensional space,
    /// find the distance to that point along the Hilbert curve.
    /// That distance will be transposed; broken into pieces and distributed into an array.
    /// This performs the most important half of the **forward Hilbert Transform**.
    /// 
    /// The number of dimensions is inferred from the length of the hilbert_axes array.
    /// 
    ///   - hilbert_axes - Point in N-space.
    ///   - bits - Depth of the Hilbert curve. If bits is one, this is the top-level Hilbert curve.
    ///   - returns The Hilbert distance (or index) as a transposed Hilbert index. 
    /// (It must be subsequently untransposed into a `BigUint` to be usable.)
    #[allow(non_snake_case)]
    pub fn hilbert_index_transposed(hilbert_axes : &[u32], bits : usize) -> Vec<u32>
    {
        let mut X = hilbert_axes.to_vec();
        let n = hilbert_axes.len(); // n: Number of dimensions
        let M : u32 = 1_u32 << (bits - 1);
        let (mut P, mut Q, mut t) : (u32, u32, u32);

        // Inverse undo
        Q = M;
        while Q > 1
        {
            P = Q - 1;

            // Split out first iteration of loop and store X[0] in a local variable to reduce array accesses.
            // Plus, we do not need to XOR X[0] with t twice since they cancel each other out.
            let mut X_0 = X[0];
            if (X_0 & Q) != 0 {
                X_0 ^= P; // invert
            }
            for i in 1..n
            {
                let X_i = X[i];
                if (X_i & Q) != 0 {
                    X_0 ^= P; // invert
                }
                else
                {
                    t = (X_0 ^ X_i) & P;
                    X_0 ^= t;
                    X[i] = X_i ^ t;
                }
            }
            X[0] = X_0;

            Q >>= 1;
        } // exchange
        // Gray encode
        let mut X_i_minus_1 = X[0];
        for i in 1..n {
            X[i] ^= X_i_minus_1;
            X_i_minus_1 = X[i];
        }
        t = 0;
        Q = M;
        while Q > 1 {
            if (X[n - 1] & Q) != 0 {
                t ^= Q - 1;
            }
            Q >>= 1;
        }
        for i in 0..n {
            X[i] ^= t;
        }
        X
    }


    /// Given the axes (coordinates) of a point in N-Dimensional space, 
    /// find the 1-dimensional distance to that point along the Hilbert curve as a `BigUint`.
    /// 
    /// The number of dimensions is the length of the `hilbert_axes` array.
    /// 
    /// Note: This is the most caller-friendly transformation from N-space to Hilbert distance.
    /// 
    ///   - `hilbertAxes` - Point in N-space.
    ///   - `bits` - Depth of the Hilbert curve. If bits is one, this is the top-level Hilbert curve.
    ///   - `interleaver_option` - Optional Interleaver that can be used to optimize the transformation. 
    ///     If only transforming a single point, this is counterproductive. If transforming many, it will speed things up.
    ///   - **Returns** - The Hilbert distance (or index) as a `BigUint`.
    pub fn hilbert_index(hilbert_axes : &[u32], bits : usize, interleaver_option : Option<&Interleaver>) -> BigUint
    {
        let transposed_index = hilbert_index_transposed(hilbert_axes, bits);
        let index = untranspose(&transposed_index, bits, interleaver_option);
        index
    }

}

#[cfg(test)]
mod tests {
    use num::BigUint;
    use std::ops::Sub;

    #[allow(unused_imports)]
    use spectral::prelude::*;
    use super::fast_hilbert;

   #[test]
    fn unpack_big_integer() {
        let big : BigUint = 329_u32.into(); // 329 = 256 + 64 + 8 + 1 = 0000101001001 (high to low) or 1001001010000 (low to high)
        let actual_bit_array = fast_hilbert::unpack_big_integer(&big, 12);
        let expected_bit_array : Vec<u8> = vec![1,0,0,1,0,0,1,0,1,0,0,0];
        asserting!("Correct number of bits").that(&actual_bit_array.len()).is_equal_to(12);
        asserting!("Correct ordering of bits").that(&actual_bit_array).is_equal_to(expected_bit_array);
    }

    #[test]
    fn uninterleave() {
        // 0110010001001100 binary = 25676 in decimal
        // .1..0..0..0..1.. => 10001 = 17
        // ..1..1..0..0..0. => 11000 = 24
        // ...0..0..1..1..0 => 00110 = 6
        // High-order bit is a pad bit, because we only need 15 bits, and two bytes hold 16
        let gray_code : BigUint = 25676_u32.into();
        let actual_axes = fast_hilbert::uninterleave(&gray_code, 5, 3, true);
        let expected_axes : Vec<u32> = vec![17, 24, 6];
        asserting("Correct uninterleave result").that(&actual_axes).is_equal_to(expected_axes);
    }

    #[test]
    fn uninterleave_2_dims_3_bits() {
        // 00001000 binary = 8 in decimal
        // ..0.1.0. => 010 = 2
        // ...0.0.0 => 000 = 0
        // 2 High-order bits are pad bits, because we only need 6 bits, and one bytes holds 8

        let gray_code : BigUint = 8_u32.into();
        let actual_axes = fast_hilbert::uninterleave(&gray_code, 3, 2, true);
        let expected_axes : Vec<u32> = vec![2, 0];
        asserting("Correct uninterleave result").that(&actual_axes).is_equal_to(expected_axes);
    }

    #[test]
    fn untranspose() {
        // 3-D Axes in decimal: [17, 24, 6]
        // 3-D Axes in Binary, 5-bits per number: [10001, 11000, 00110]
        // BigUint, 15 bit + 1 pad binary number: 0110010001001100 = 25676 in decimal
        let axes : Vec<u32> = vec![17, 24, 6];
        let actual = fast_hilbert::untranspose(&axes, 5, None);
        let expected : BigUint = 25676_u32.into();
        asserting("Correct untranspose result").that(&actual).is_equal_to(expected);
    }

    #[test]
    fn interleave_be() {
        // 3-D Axes in Binary, 5-bits per number: decimal [17, 24, 6] = binary [10001, 11000, 00110]
        // Distributed to two bytes, putting one bit padding in the high-bit position:
        //    [01100100,01001100] = [100,76] in decimal
        let axes : Vec<u32> = vec![17, 24, 6];
        let actual = fast_hilbert::interleave_be(&axes, 5, None);
        let expected : Vec<u8> = vec![100,76];
        asserting("Correct interleave result").that(&actual).is_equal_to(expected);
    }

    fn abs_difference<T: Sub<Output = T> + Ord>(x: T, y: T) -> T {
        if x < y {
            y - x
        } else {
            x - y
        }
    }

    /// Compare two vectors to see if they differ in only one dimension and that difference is exactly one.
    fn verify_difference(vec1 : &[u32], vec2 : &[u32]) -> Result<usize, String> {
        let mut changed_dimension = None;
        for dimension in 0..(vec1.len()) {
            let diff = abs_difference(vec1[dimension], vec2[dimension]);
            if diff != 0 {
                if changed_dimension.is_some() {
                    return Err(format!("More than one dimension changed: {} and {}. Vectors: {:?} and {:?}", dimension, changed_dimension.unwrap(), vec1, vec2));
                }
                else {
                    changed_dimension = Some(dimension);
                    if diff != 1 {
                        return Err(format!("Dimension {} changed by {} instead of 1", dimension, diff));
                    }
                }
            }
        }
        match changed_dimension {
            Some(dim) => Ok(dim),
            None => Err("No dimensions changed".to_string())
        }
    }

    /// Verifying that the Hilbert transformation is correct in high dimensions is impractical. 
    /// This test verifies three properties of the transformation for a curve with 5 bits per coordinate in 3 dimensions: 
    /// 
    ///   1. It is **invertible**. Performing the transform of a Hilbert index to a vector then the reverse transform should
    ///      yield the original index.
    ///   2. It is **exhaustive**. Every possible Hilbert Index is visited, from zero to 2^(Bits*Dimensions) - 1
    ///   3. It takes **baby steps**. Each consecutive vector yielded by increasing the Hilbert Index by one 
    ///      should differ from the previous vector in a single dimension, and the change in coordinate value 
    ///      should only be +1 or -1.
    ///   4. It is **twisty**. The dimension that changes in one iteration should _seldom_ match the dimension 
    ///      that changed in the immediately previous iteration. In rare instances, I have seen it move the same 
    ///      direction three times in a row. I do not know if it ever moves farther than that in the same direction
    ///      without turning. 
    /// 
    /// This tests 32,768 forward and backwards transforms (for index = 0 to index = 32767, inclusive). 
    #[test]
    fn hilbert_round_trip(){
        let bits = 5;
        let dimensions = 3;
        let mut previous : Option<(Vec<u32>, Option<usize>)> = None;
        let mut same_dimension_changed_twice = false;
        let mut same_dimension_changed_thrice = false;

        // Iterating over all possible Hilbert index values without a panic plus passing the invertibility test 
        // satisfies the "exhaustive" requirement.
        for i in 0..32768_u32 {
            let expected_hilbert_index : BigUint = i.into();
            let coordinates = fast_hilbert::hilbert_axes(&expected_hilbert_index, bits, dimensions);
            let actual_hilbert_index = fast_hilbert::hilbert_index(&coordinates, bits, None);

            asserting(&format!("Invertible for i = {}", i)).that(&actual_hilbert_index).is_equal_to(expected_hilbert_index);
            previous = match previous {
                None => Some((coordinates.clone(), None)),
                Some((previous_coordinates, None)) => {
                    // verify_difference ensures that only one dimension changed and the coordinate change was +-1. 
                    same_dimension_changed_twice = false;
                    same_dimension_changed_thrice = false;
                    match verify_difference(&coordinates, &previous_coordinates) {
                        Ok(changed_dim) => {
                            Some((coordinates, Some(changed_dim)))
                        },
                        Err(msg) => {
                            panic!(format!("Error: {} for Hilbert index = {}", msg, i));
                        }
                    }
                },
                Some((previous_coordinates, Some(previous_changed_dim))) => {
                    // verify_difference ensures that only one dimension changed and the coordinate change was +-1. 
                    match verify_difference(&coordinates, &previous_coordinates) {
                        Ok(changed_dim) => {
                            if previous_changed_dim == changed_dim {
                                asserting(&format!("Dimension {} changed more than thrice in a row at i = {}", previous_changed_dim, i))
                                    .that(&same_dimension_changed_thrice).is_equal_to(false);
                                if same_dimension_changed_twice {
                                    same_dimension_changed_thrice = true;
                                }
                                else {
                                    same_dimension_changed_twice = true;
                                }
                            }
                            else {
                                same_dimension_changed_twice = false;
                                same_dimension_changed_thrice = false;
                            }
                            Some((coordinates, Some(changed_dim)))
                        },
                        Err(msg) => {
                            asserting(&format!("Error: {} for Hilbert index = {}", msg, i)).that(&true).is_equal_to(false);
                            panic!("Should never reach this line");
                        }
                    }
                }
            }
        }
    }

    /// Verify that a simple case of two dimensions and 3 bits per dimension yields the exact correct solution.
    /// 
    /// There is not only one way to arrange points in a Hilbert curve, since they can be
    /// rotations and reflections of one another. This test assumes a certain solution,
    /// derived from the C# implementation.
    #[test]
    fn exact_2_dimensions_3_bits() {
        let bits : usize = 3;
        let dimensions : usize = 2;
        let num_points : usize = 1 << (bits * dimensions);
        let expected : Vec<Vec<u32>> = vec![
            vec![0,0], vec![0,1], vec![1,1], vec![1,0], vec![2,0], vec![3,0],
            vec![3,1], vec![2,1], vec![2,2], vec![3,2], vec![3,3], vec![2,3],
            vec![1,3], vec![1,2], vec![0,2], vec![0,3], vec![0,4], vec![1,4],
            vec![1,5], vec![0,5], vec![0,6], vec![0,7], vec![1,7], vec![1,6],
            vec![2,6], vec![2,7], vec![3,7], vec![3,6], vec![3,5], vec![2,5],
            vec![2,4], vec![3,4], vec![4,4], vec![5,4], vec![5,5], vec![4,5],
            vec![4,6], vec![4,7], vec![5,7], vec![5,6], vec![6,6], vec![6,7],
            vec![7,7], vec![7,6], vec![7,5], vec![6,5], vec![6,4], vec![7,4],
            vec![7,3], vec![7,2], vec![6,2], vec![6,3], vec![5,3], vec![4,3],
            vec![4,2], vec![5,2], vec![5,1], vec![4,1], vec![4,0], vec![5,0],
            vec![6,0], vec![6,1], vec![7,1], vec![7,0]
        ];
        for i in 0..num_points {
            let hilbert_index : BigUint = i.into();
            let coordinates = fast_hilbert::hilbert_axes(&hilbert_index, bits, dimensions);
            asserting(&format!("Hilbert Index = {}. Expected {:?}. Actual {:?}", i, expected[i], coordinates)).that(&coordinates).is_equal_to(expected[i].clone());
        }
    }

    #[test]
    fn single_transform_dim2_bits3_index8() {
        let bits : usize = 3;
        let dimensions : usize = 2;
        let index = 8_u32;
        let hilbert_index : BigUint = index.into();
        let actual_point = fast_hilbert::hilbert_axes(&hilbert_index, bits, dimensions);
        let expected_point = vec![2,2];
        asserting(&format!("Hilbert Index = {}. Expected {:?}. Actual {:?}", index, expected_point, actual_point)).that(&actual_point).is_equal_to(expected_point);
    }

    
    fn dump_hilbert(dimensions : usize, bits : usize, from : u32, to : u32){
        println!("{} dimensions, {} bits, i = {} to {}", dimensions, bits, from, to);
        for i in from..to {
            let expected_hilbert_index : BigUint = i.into();
            let coordinates = fast_hilbert::hilbert_axes(&expected_hilbert_index, bits, dimensions);
            println!("{} => {:?}", i, coordinates);
        }

    }

    #[test]
    #[ignore]
    fn show_hilbert(){
        dump_hilbert(3, 5, 0, 40);
        panic!("show_hilbert");
    }
}
