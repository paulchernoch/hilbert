//! The `point_list` module can convert point-like data in forms defined by the caller
//! into `Point` structs. It defines the functions `make_points_i32` and `make_points_f64`. 
//! They can detect the range of values in the input data and the number of bits required
//! to code those values without losing precision, then translate and scale those values
//! in the process of creating `Points`. `Points` created in this way are guaranteed
//! to be acceptable to the **Hilbert Curve** transformation.
use super::point::Point;
use super::normalize::{ IntegerDataRange, FloatDataRange };

/// Take unnormalized `i32` data, normalize it, and create Point objects 
/// (which hold unsigned values). 
/// 
///    - `points` - The point data to be normalized and optionally scaled. 
///    - `starting_id` - The first point created will get this value for id and each one after that will increment by one.
///    - `range_option` - If supplied, range is used in the normalization, 
///                       otherwise a range is calculated by studying all the input points. 
///    - `bits_allocated` - If supplied, compress the values to this number of bits, 
///                         otherwise use the minimum number of bits required to faithfully represent the full range of data. 
///    - **return** - A tuple holding a Vec of Points whose coordinate values have been normalized and optionally scaled, 
///                   and the number of bits used to represent each dimension. (The number of bits used must be fed to the Hilbert transformation.)
pub fn make_points_i32<I>(points : &[I], starting_id : usize, range_option : Option<IntegerDataRange>, bits_allocated : Option<usize>) -> (Vec<Point>, usize)
// The "for<'a>" below is "higher-ranked lifetime" notation. Got it from stackoverflow. Confused about what it does but it is needed!
where
    for<'a> &'a I: IntoIterator<Item = &'a i32>, 
{
    let mut next_id = starting_id;
    let mut converted_points = Vec::with_capacity(points.len());
    let mut coordinates = Vec::new();

    // If a range was not passed in, we must infer it by studying the data. 
    let range = match range_option {
        Some(given_range) => given_range.clone(),
        None => IntegerDataRange::from_i32(points)
    };
    let bits_used = match bits_allocated {
        Some(bits) => bits,
        None => range.bits_required
    };
    for input_point in points {
        coordinates.clear();
        for coordinate in input_point {
            let transformed_coordinate = match bits_allocated {
                Some(bits) => range.compress(*coordinate, bits),
                None => range.normalize(*coordinate)
            };
            coordinates.push(transformed_coordinate);
        }
        converted_points.push(Point::new(next_id, &coordinates));
        next_id += 1;
    }
    (converted_points, bits_used)
}

/// Take unnormalized `f64` data, normalize it, and create Point objects 
/// (which hold unsigned integer values). 
/// 
///    - `points` - The point data to be normalized and optionally scaled. 
///    - `starting_id` - The first point created will get this value for id and each one after that will increment by one.
///    - `range_option` - If supplied, range is used in the normalization, 
///                       otherwise a range is calculated by studying all the input points. 
///    - `bits_allocated` - If supplied, compress the values to this number of bits, 
///                         otherwise use the minimum number of bits required to faithfully represent the full range of data. 
///    - `scale` - Multiply all coordinates by this before rounding to integers. 
///    - **return** - A tuple holding a Vec of Points whose coordinate values have been normalized and optionally scaled, 
///                   and the number of bits used to represent each dimension. (The number of bots used must be fed to the Hilbert transformation.)
pub fn make_points_f64<I>(points : &[I], starting_id : usize, range_option : Option<FloatDataRange>, bits_allocated : Option<usize>, scale : f64) -> (Vec<Point>, usize)
// The "for<'a>" below is "higher-ranked lifetime" notation. Got it from stackoverflow. Confused about what it does but it is needed!
where
    for<'a> &'a I: IntoIterator<Item = &'a f64>, 
{
    let mut next_id = starting_id;
    let mut converted_points = Vec::with_capacity(points.len());
    let mut coordinates = Vec::new();

    // If a range was not passed in, we must infer it by studying the data. 
    let range = match range_option {
        Some(given_range) => given_range.clone(),
        None => FloatDataRange::from_f64(points, scale)
    };
    let bits_used = match bits_allocated {
        Some(bits) => bits,
        None => range.bits_required
    };
    for input_point in points {
        coordinates.clear();
        for coordinate in input_point {
            let transformed_coordinate = match bits_allocated {
                Some(bits) => range.compress(*coordinate, bits),
                None => range.normalize(*coordinate)
            };
            coordinates.push(transformed_coordinate);
        }
        converted_points.push(Point::new(next_id, &coordinates));
        next_id += 1;
    }
    (converted_points, bits_used)
}

#[cfg(test)]
mod tests {

    #[allow(unused_imports)]
    use spectral::prelude::*;
    use crate::point::Point;

    /// Make some points and verify that they are properly normalized and the required number of bits is properly identified.
    #[test]
    fn make_points_f64() {
        let float_points = test_points_f64();
        let (actual_points, bits_used) = super::make_points_f64(&float_points, 0, None, None, 10.0);
        asserting("correct bits_used").that(&bits_used).is_equal_to(11);
        let expected_points = vec![
            Point::new(0, &vec![895, 1053, 1037]),
            Point::new(0, &vec![952, 1202, 2002]),
            Point::new(0, &vec![1420, 0, 1000])
        ];
        for i in 0..3 {
            let expected_point = &expected_points[i];
            let actual_point = &actual_points[i];
            asserting(&format!("For point {:?}", expected_point)).that(&expected_point.are_coordinates_identical(actual_point)).is_equal_to(true);
        }
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
