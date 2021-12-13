//! The `normalize` module contains `IntegerDataRange` and `FloatDataRange`.
//! These structs are used to normalize point coordinates to prepare them for the **Hilbert transform**.

use std::{i32, u32, u64, f64};
use std::iter::{Iterator, IntoIterator};




/// Smallest power of two such that two raised to that power is greater than or equal to the given number.
///
///   - `n` - Find the smallest power of two that equals or exceeds this value.
///   - Returns - The power, not the number (unlike u64::next_power_of_two, which returns the number, not its log base 2). 
/// 
/// Examples: 
///   - smallest_power_of_two(7) returns 3, because 2^3 = 8 which is **greater than** 7.
///   - smallest_power_of_two(16) returns 4, because 2^4 = 16 which is **equal to** 16.
pub fn smallest_power_of_two<N>(n : N) -> usize
where N : Into<u64>
{
    let n64 : u64 = n.into();
    let mut r = 1_u64;
    let mut log_two = 0;
    while r < n64
    {
        r <<= 1;
        log_two += 1;
    }
    log_two
}

/// Find the minimum number of bits required to represent the given number.
/// 
/// Examples: 
///   - bits_required(7) returns 3, because 2^3 = 8 which is **greater than** 7.
///   - bits_required(16) returns 5, because 2^5 = 32 which is **greater** 16.
pub fn bits_required<N>(n : N) -> usize
where N : Into<u64>
{
    smallest_power_of_two(n.into() + 1)
}


/// Holds the observed range of values seen in a set of points whose coordinates are 32-bit integers (signed or unsigned).
/// 
/// From this information we can determine how to translate the coordinates so that their normalized range has a low of zero.
/// This permits us to use the fewest number of bits to represent them when constructing a **Hilbert index**.
/// 
///   - Use the `normalize` method to adjust a coordinate value into range such that the full precision of the number is preserved.
///   - Use the `compress` method to both normalize and shrink the coordinate value, 
///     to reduce the number of bits used per coordinate, at the expense of less precision.
/// 
/// **Motivation**. 
/// 
///   1. The Hilbert Curve transformation requires _non-negative numbers_, so all values must be translated until 
///      the lowest is zero. That means that if the lowest value is negative, we shift to the positive direction, 
///      or contrariwise to the negative direction. 
///   2. The _execution time_ of the Hilbert transformation is directly proportional to the _number
///      of bits_ used to represent each value. Thus if we can sacrifice some precision, each value can be
///      right-shifted by the same number of bits if we wish to compress the range. 
#[derive(Clone, PartialEq, Debug)]
pub struct IntegerDataRange {
    /// Lowest value of any coordinate of any point in a collection of points.
    pub low : i64,

    /// Highest value of any coordinate of any point in a collection of points.
    pub high : i64,

    /// Minimum number of bits required to represent a normalized value without loss of information. 
    /// 
    /// Examples: 
    /// 
    ///    - If `low` is -10 and `high` is 24 then the range is 34 so 6 bits are required to represent all values in that range. 
    ///    - If `low` is 50 and `high` is 67 then the range is 17 so 5 bits are required to represent all values in that range. 
    pub bits_required : usize
}

impl IntegerDataRange {

    /// Create an `IntegerDataRange` without reference to particular data. 
    pub fn new<I>(low_i : I, high_i : I) -> Self 
    where I : Into<i64>
    {
        let low = low_i.into();
        let high = high_i.into();
        IntegerDataRange { 
            low, 
            high,
            bits_required : bits_required((high - low) as u64)
        }
    }

    /// Study all `i32` coordinates in all points to find the minimum and maximum values. 
    pub fn from_i32<I>(points : &[I]) -> Self 
    // The "for<'a>" below is "higher-ranked lifetime" notation. Got it from stackoverflow. Confused about what it does but it is needed!
    where
        for<'a> &'a I: IntoIterator<Item = &'a i32>, 
    {
        let mut low = i32::MAX;
        let mut high = i32::MIN;
        for point in points.iter() {
            for coordinate in point.into_iter().map(|c| *c) {
                if coordinate > high { high = coordinate; }
                if coordinate < low { low = coordinate; }
            }
        }
        Self::new(low, high)
    }

    /// Study all `u32` coordinates in all points to find the minimum and maximum values. 
    pub fn from_u32<I>(points : &[I]) -> Self 
    where
        for<'a> &'a I: IntoIterator<Item = &'a u32>,
    {
        let mut low = u32::MAX as i64;
        let mut high = 0_i64;
        for point in points.iter() {
            for coordinate in point.into_iter().map(|c| *c as i64) {
                if coordinate > high { high = coordinate; }
                if coordinate < low { low = coordinate; }
            }
        }
        Self::new(low, high)
    }

    /// Range from low to high value.
    pub fn range(&self) -> u32 { (self.high - self.low) as u32 }

    /// Normalize an integer convertible to an i64, shifting it enough so that the minimum value found in any point
    /// is shifted to zero and the maximum value is shifted to `range`, and using the full number of bits required for the range.
    pub fn normalize<I>(&self, coordinate : I) -> u32 
    where I : Into<i64>
    {
        (coordinate.into() - self.low) as u32
    }

    /// Normalize a `coordinate` value, shifting it enough so that the minimum value found in any point
    /// is shifted to zero and the maximum value is shifted to `range`, then optionally compressing the range by bit shifting
    /// such that no more than the given number of bits are required for the largest value.
    pub fn compress<I>(&self, coordinate : I, bits_allocated : usize) -> u32 
    where I : Into<i64>
    {
        let uncompressed = self.normalize(coordinate);
        if bits_allocated < self.bits_required { uncompressed >> (self.bits_required - bits_allocated) }
        else { uncompressed }
    }

}

/// The observed range of values seen in a set of points whose coordinates are signed 64-bit floats.
/// 
/// From this information we can determine how to translate the coordinates so that their normalized range has a low of zero.
/// This permits us to use the fewest number of bits to represent them when constructing a **Hilbert index**.
/// 
/// **Motivation**. 
/// 
///   1. The Hilbert Curve transformation requires _non-negative numbers_, so all values must be translated until 
///      the lowest is zero. That means that if the lowest value is negative, we shift to the positive direction, 
///      or contrariwise to the negative direction. 
///   2. The _execution time_ of the Hilbert transformation is directly proportional to the _number
///      of bits_ used to represent each value. Thus if we can sacrifice some precision, each value can be
///      multiplied by a scale factor and rounded to the nearest integer to compress the value (and sacrifice information). 
#[derive(Clone, PartialEq, Debug)]
pub struct FloatDataRange {
    /// Lowest value of any coordinate of any point in a collection of points.
    pub low : f64,


    /// Highest value of any coordinate of any point in a collection of points.
    pub high : f64,

    /// Multiplier to apply before rounding to an integer value, sacrificing some precision. 
    /// For example, if you want to encode values such that you preserve precision to the hundredth place, `scale` should be 100.
    pub scale : f64, 

    /// Minimum number of bits required to represent a normalized value without loss of information for the given scale factor. 
    /// 
    /// Example: 
    /// 
    ///    - If `low` is -5.02 and `high` is 3.13 then the range is 8.15. If `scale` is 100, the range becomes 815, 
    /// so 10 bits are required to represent all values in that range. 
    pub bits_required : usize
}

impl FloatDataRange {

    /// Create a `FloatDataRange` without reference to particular data. 
    pub fn new(low : f64, high : f64, scale : f64) -> Self {
        let scaled_range = ((high - low) * scale).ceil() as u32;
        FloatDataRange { 
            low, 
            high,
            scale,
            bits_required : bits_required(scaled_range)
        }
    }

    /// Study all `f64` coordinates in all points to find the minimum and maximum values. 
    pub fn from_f64<I>(points : &[I], scale : f64) -> Self 
    // The "for<'a>" below is "higher-ranked lifetime" notation. Got it from stackoverflow. Confused about what it does but it is needed!
    where
        for<'a> &'a I: IntoIterator<Item = &'a f64>, 
    {
        let mut low = f64::MAX;
        let mut high = f64::MIN;
        for point in points.iter() {
            for coordinate in point.into_iter().map(|c| *c) {
                if coordinate > high { high = coordinate; }
                if coordinate < low { low = coordinate; }
            }
        }
        FloatDataRange::new(low, high, scale)
    }

    /// Range from low to high value.
    pub fn range(&self) -> f64 { self.high - self.low }

    /// Normalize an `f64` `coordinate` value, shifting it enough so that the minimum value found in any point
    /// is shifted to zero and the maximum value is shifted to `range`, and using the full number of bits required for the range
    /// multiplied by the scale.
    pub fn normalize(&self, coordinate : f64) -> u32 {
        ((coordinate - self.low) * self.scale).ceil() as u32
    }

    /// Normalize an `f64` `coordinate` value, shifting it enough so that the minimum value found in any point
    /// is shifted to zero and the maximum value is shifted to `range`, then optionally compressing the range by bit shifting
    /// such that no more than the given number of bits are required for the largest value.
    pub fn compress(&self, coordinate : f64, bits_allocated : usize) -> u32 {
        let uncompressed = self.normalize(coordinate);
        if bits_allocated < self.bits_required { uncompressed >> (self.bits_required - bits_allocated) }
        else { uncompressed }
    }

}

#[cfg(test)]
mod tests {

    #[allow(unused_imports)]
    use spectral::prelude::*;
    use super::{IntegerDataRange, FloatDataRange};

    #[test]
    fn from_i32() {
        let points = test_points_i32();
        let actual_range = IntegerDataRange::from_i32(&points);
        let expected_range = IntegerDataRange::new(-100, 100);
        asserting("Range correct").that(&actual_range).is_equal_to(expected_range);
    }

    #[test]
    fn normalize_i32() {
        let points = test_points_i32();
        let range = IntegerDataRange::from_i32(&points);
        let actual_normalized = range.normalize(-16);
        let expected_normalized = 84;
        asserting("Normalized value correct").that(&actual_normalized).is_equal_to(expected_normalized);
    }    

    #[test]
    fn compress_i32() {
        let points = test_points_i32();
        let range = IntegerDataRange::from_i32(&points);
        let actual_compressed = range.compress(-16, 7);
        let expected_compressed = 42;
        asserting("Compressed value correct").that(&actual_compressed).is_equal_to(expected_compressed);
    }

    #[test]
    fn smallest_power_of_two() {
        let in_out_pairs = vec![(0_u64,0), (1,0), (2,1), (3,2), (4,2), (5,3), (6,3), (7,3), (8,3), (9,4)];
        for (n, expected) in in_out_pairs {
            let actual = super::smallest_power_of_two(n);
            asserting(&format!("n = {}, actual = {}, expected = {}", n, actual, expected)).that(&actual).is_equal_to(expected);
        }
    }

    #[test]
    fn from_f64() {
        let points = test_points_f64();
        let actual_range = FloatDataRange::from_f64(&points, 10.0);
        let expected_range = FloatDataRange::new(-100.0, 100.19, 10.0);
        asserting("Range correct").that(&actual_range).is_equal_to(expected_range);
    }

    #[test]
    fn normalize_f64() {
        let points = test_points_f64();
        let range = FloatDataRange::from_f64(&points, 10.0);
        let actual_normalized = range.normalize(-16.0);
        let expected_normalized = 840;
        asserting("Normalized value correct").that(&actual_normalized).is_equal_to(expected_normalized);
    }

    /// Compress a value in a range that requires 11 bits down to 7 bits, forcing a division by 16 and loss of precision. 
    #[test]
    fn compress_f64() {
        let points = test_points_f64();
        let range = FloatDataRange::from_f64(&points, 10.0);
        let actual_compressed = range.compress(-16.0, 7);
        let expected_compressed = 52;
        asserting("Compressed value correct").that(&actual_compressed).is_equal_to(expected_compressed);
    }    

    fn test_points_i32() -> Vec<Vec<i32>> {
       vec![
           vec![-10, 5, 3],
           vec![-4, 20, 100],
           vec![42, -100, 0]
       ]
    }

    fn test_points_f64() -> Vec<Vec<f64>> {
       // Range is 100.9 - -100.0 => 200.9
       // If scaled by 10, the range is 2009. 
       vec![
           vec![-10.5, 5.27, 3.66],
           vec![-4.802, 20.2, 100.19],
           vec![42.0, -100.0, 0.0]
       ]
    }

}