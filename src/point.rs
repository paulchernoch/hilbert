//! The `point` module defines the `Point` struct, which represents an N-dimensional point.
//! These `Points` are suitable for the Hilbert Transformation, because they do not have negative values
//! and the total number of bits required to represent each coordinate is controlled by the construction process.
//! 
//! For the most convenient means of constructing many `Points`, see the `point_list` module.
use std::hash::{Hash, Hasher};
use num::BigUint;
use super::transform::fast_hilbert;
use super::permutation::Permutation;

/// An immutable, N-dimensional point with unsigned integer coordinates and an optimized distance function.
/// 
/// If the `Point` was created using `make_points_i32` or `make_points_f64` from the `point_list` module, 
/// or if the caller took special care to ensure the coordinate values are in the desired range for the **Hilbert Curve**,
/// then the `Point` is fit to be converted into a **Hilbert Index** using the `hilbert_transform` method.
/// 
/// Use `hilbert_sort` if you are uninterested in the **Hilbert Index** and merely need many points sorted by the index.
/// 
/// # Examples
///
/// ```
///         use crate::hilbert::point::Point;
///         use crate::hilbert::point_list;
/// 
///         // 1. Create two 3-D points and get the square of the distance between them.
///         let p1 = Point::new(0, &[3, 4, 5]);
///         let p2 = Point::new(1, &[0, 8, 10]);
///         let sqr_dist = p1.square_distance(&p2);
///         assert!(sqr_dist == 50, "Square distance should be 50");
/// 
///         // 2. Perform the Hilbert Transform on a single point,
///         //    using 5 bits per dimension (which assumes no coordinate exceeds 31).
///         let index1 = p1.hilbert_transform(5);
/// 
///         // 3. Create several points and normalize them.
///         //    This will ensure that the ids begin at zero and that all values
///         //    are multiplied by 10.0 before being rounded to the nearest integer,
///         //    to preserve the predetermined precision.
///         let point_data : Vec<Vec<f64>> = vec![
///            vec![-10.5, 5.27, 3.66],
///            vec![-4.802, 20.2, 100.19],
///            vec![42.0, -100.0, 0.0]
///         ];
///         let (mut points, bits) = point_list::make_points_f64(&point_data, 0, None, None, 10.0);
/// 
///         // 4. Sort the points by the Hilbert Curve, using 11 bits per dimension,
///         //    because the range of data is 200.19, multiplied by the scale of 10 
///         //  yields 2001.9, ceiling of that yields 2002, which is between 1024 (2^10) 
///         //  and 2048 (2^11), so 11 bits are required to store the 
///         //  highest coordinate value.
///         Point::hilbert_sort(&mut points, 11);
/// ```
#[derive(Debug, Clone)]
pub struct Point {
    /// Presumably unique id.
    id : usize,

    /// Coordinate values for point.
    coordinates : Vec<u32>,

    /// Square of the distance from the origin to the `Point`.
    square_magnitude : u64
}

impl Point {
    /// Construct a new Point with the given id. 
    /// 
    /// Note: It is the caller's duty to ensure that ids remain unique. 
    pub fn new(id : usize, coords : &[u32]) -> Point {
        Point {
            id, 
            coordinates : coords.iter().map(|c| c.clone()).collect(), 
            square_magnitude : coords.iter().fold(0_u64, |sum,coord| sum + (coord.clone() as u64) * (coord.clone() as u64))
        }
    }

    /// Apply a permutation to the `Point` to generate a new `Point` 
    /// where the coordinate values have been reordered, but the `id` is the same. 
    pub fn permute(&self, permutation : &Permutation) -> Self {
        Point::new(self.id, &permutation.apply(&self.coordinates))
    }

    /// Number of dimensions for the Point. 
    pub fn dimensions(&self) -> usize { 
        self.coordinates.len() 
    }

    /// Gets the `id` for the Point. 
    pub fn get_id(&self) -> usize { self.id }

    /// Gets the coordinate values for the point.
    pub fn get_coordinates(&self) -> &Vec<u32> { &self.coordinates }

    /// Square of the Cartesian distance between the two points. 
    /// This is optimized by using the following algebra: 
    ///
    /// ```ignore
    /// 
    ///  //    (A - B)² = A² + B² - 2·A·B
    /// ```
    /// If you extend this to entire vectors, then:
    ///    - A² is the square magnitude of vector A, 
    ///    - B² is the square magnitude of vector B,
    ///    - 2·A·B is twice the dot product of A and B. 
    /// This saves about N subtractions.  
    pub fn square_distance(&self, other_point : &Point) -> u64 {
        if self.dimensions() != other_point.dimensions() {panic!()}
        let dotprod = self
            .coordinates
            .iter()
            .zip(other_point.coordinates.iter())
            .map(|(&a, &b)| (a as u64)*(b as u64))
            .sum();
        self.square_magnitude + other_point.square_magnitude - 2*dotprod
    }

    /// Compares the coordinates of the two `Points` to see if they match in value and sequence. 
    /// 
    /// This comparison ignores the ids of the `Points`. 
    pub fn are_coordinates_identical(&self, other_point : &Point) -> bool {
        self.square_magnitude == other_point.square_magnitude && self.coordinates == other_point.coordinates
    }

    /// Perform a Hilbert transform from a single N-dimensional point to a 1-Dimensional index. 
    /// 
    ///   - `bits_per_dimension` - Number of bits used to encode each coordinate of the Point. 
    ///      This value must be at least one.
    ///      The `BigUint` Hilbert Index for the `Point` will be composed of enough bytes
    ///      to hold N * `bits_per_dimension` bits, where N is the number of points. 
    ///      If any coordinate of the `Point` has a value c >= 2^ `bits_per_dimension`, 
    ///      results are undefined and will be incorrect.
    /// 
    /// Note: This should not be called for `Points` with fewer than two dimensions. 
    pub fn hilbert_transform(&self, bits_per_dimension : usize) -> BigUint {
        fast_hilbert::hilbert_index(&self.coordinates, bits_per_dimension, None)
    }

    /// Sort a collection of `Points` in ascending **Hilbert Index** order. 
    /// 
    /// This method discards the Hilbert indices when done.
    ///    
    ///   - `points` - Points to sort. 
    ///   - `bits_per_dimension` - Number of bits used to encode each dimension. 
    ///      The `BigUint` Hilbert Index for each point will be composed of enough bytes
    ///      to hold N * `bits_per_dimension` bits, where N is the number of points. 
    ///      If any coordinate of any `Point` has a value c >= 2^ `bits_per_dimension`, 
    ///      results are undefined and will be incorrect.
    /// 
    /// The points are sorted in place.
    pub fn hilbert_sort(points : &mut Vec<Point>, bits_per_dimension : usize) {
        points.sort_by_cached_key(|point| point.hilbert_transform(bits_per_dimension));
    }

    /// Sort a collection of `Points` in ascending **Hilbert Index** order after applying a
    /// `Permutation` to each point. This has the effect of generating a rotated or reflected
    /// Hilbert curve. For N-dimensional points, there are N! variations of the Hilbert curve that can
    /// be generated in this way. One curve among these may be more useful than another for a given purpose.
    /// 
    /// This method discards the Hilbert indices when done.
    ///    
    ///   - `points` - Points to sort. 
    ///   - `bits_per_dimension` - Number of bits used to encode each dimension. 
    ///      The `BigUint` Hilbert Index for each point will be composed of enough bytes
    ///      to hold N * `bits_per_dimension` bits, where N is the number of points. 
    ///      If any coordinate of any `Point` has a value c >= 2^ `bits_per_dimension`, 
    ///      results are undefined and will be incorrect.
    ///   - `permutation` - Defines how to reorder the many coordinates in a repeatable way.
    pub fn hilbert_sort_permuted(points : &mut Vec<Point>, bits_per_dimension : usize, permutation : &Permutation) {
        points.sort_by_cached_key(|point| point.permute(permutation).hilbert_transform(bits_per_dimension));
    }
}

impl Hash for Point {
    /// To remain consistent with equals, this incorporates the `id` and `square_magnitude` into the hash. 
    fn hash<H: Hasher>(&self, state: &mut H) { 
        self.id.hash(state); 
        self.square_magnitude.hash(state);
    }
}

impl PartialEq for Point {
    /// Define equality to only compare the id and square_magnitude, for speed.
    /// 
    /// Note: To compare all the coordinate values, use `are_coordinates_identical` instead.
    fn eq(&self, other: &Self) -> bool { self.id == other.id && self.square_magnitude == other.square_magnitude }
}

#[cfg(test)]
/// Tests of the Point methods.
mod tests {
    #[allow(unused_imports)]
    use spectral::prelude::*;
    extern crate rand;
    use rand::thread_rng;
    use rand::seq::SliceRandom;
    use super::Point;
    use crate::point_list;

    #[test]
    fn square_distance() {
        let p1 = Point::new(0, &[3, 4, 5]);
        let p2 = Point::new(0, &[0, 8, 10]);
        asserting!("Square distance between points").that(&p1.square_distance(&p2)).is_equal_to(50);
    }

    #[test]
    fn hilbert_sort() {
        let expected_points = hilbert_ordered_2_dimensions_3_bits();
        let mut actual_points = expected_points.clone();
        let mut rng = thread_rng();

        // After shuffling, the points should not be in the same order.
        actual_points.shuffle(&mut rng);
        asserting("Points not sorted by Hilbert Index").that(&actual_points).is_not_equal_to(expected_points.clone());

        // After sorting, the points should once again be in teh same order.
        Point::hilbert_sort(&mut actual_points, 3);
        asserting("Points sorted by Hilbert Index").that(&actual_points).is_equal_to(expected_points);
    }

    /// Construct all sixty-four possible 2-dimensional Points, each with coordinates requiring no more than 3 bits to represent
    /// and return them in a Vec manually presorted in Hilbert Index order.
    fn hilbert_ordered_2_dimensions_3_bits() -> Vec<Point> {
        let point_data : Vec<Vec<i32>> = vec![
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
        point_list::make_points_i32(&point_data, 0, None, None).0
    }
}
