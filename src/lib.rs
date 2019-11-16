//! The Hilbert module contains these sub-modules: 
//! 
//!   - normalize - To prepare data for transformation.
//!   - permutation - To reorder point coordinates as a means of generating alternate Hilbert curves.
//!   - point_list - Prepare points from i32 or f64 data, deriving and applying a consistent transform to each.
//!   - point - Represents an N-dimensional point suitable to be transformed and sorted by the Hilbert Curve.
//!   - transform - Performs the core Hilbert Curve transform and its inverse.
extern crate num;
pub mod point;
pub mod point_list;
pub mod transform;
pub mod interleaver;
pub mod normalize;
pub mod permutation;


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
