use num::BigUint;
use crate::fast_hilbert::{hilbert_axes, hilbert_index };

/// Supplies everything needed to evaluate the integrand function for a given Hilbert index. 
/// 
///  - How to map a Hilbert Index to a point in problem coordinate space (via Hilbert Axes)
///  - How to map a point in the problem coordinate space to a Hilbert Index (via Hilbert Axes)
///  - How to evaluate integrand function on a point in problem coordinate space
///  - How to compute (possibly weighted) distance between two points in problem coordinate space
pub trait HilbertModel {
    /// Number of dimensions for points in the space.
    fn dimensions(&self) -> u16;

    /// Number of bits of precision to use for each dimension.
    fn bits(&self) -> u8;

    /// Map a point in the problem coordinate space to a Hilbert Index (via Hilbert Axes)
    fn to_hilbert_index(&self, point: &[f64]) -> num::BigUint;

    /// Map a Hilbert Index to a point in problem coordinate space (via Hilbert Axes)
    fn to_point(&self, hilbert_index: &BigUint) -> Vec<f64>;

    /// Evaluate the integrand function, returning NAN if the function is not defined for the given point.
    fn evaluate(&self, point: &[f64]) -> f64;

    /// Compute the distance between two points in the problem coordinate system, 
    /// which may be a weighted distsnce.
    fn distance(&self, point1: &[f64], point2: &[f64]) -> f64;

    /// If the class implementing this trait cached any values, free them.
    fn clear(&self);

}

/// Defines the range over which to integrate for a dimension
/// and how to scale this dimension relative to others 
/// when applying a distance metric.
#[derive(Clone,Copy,Debug)]
pub struct DimensionRange {
    /// Minimum value in the problem space for a coordinate in this dimension.
    continuous_minimum: f64,

    /// Maximum value in the problem space for a coordinate in this dimension.
    continuous_maximum: f64,

    /// One greater than the maximum permitted value in the quantized integer-valued dimension.
    quantized_maximum: u64,

    /// Derived scale factor to multiply by a problem space value to get an integer-valued coordinate.
    scale: f64,

    /// Relative weight compared to other dimensions,
    /// to be used for a weighted Euclidean distance. 
    weight: f64
}

impl DimensionRange {
    /// Construct a DimensionRange.
    pub fn new(min: f64, max: f64, weight: f64, bits: usize) -> DimensionRange {
        if max <= min {
            panic!("Maximum must exceed minimum for dimension");
        }
        if weight <= 0.0 {
            panic!("Weight must be positive");
        }
        if bits < 1 || bits > 32 {
            panic!("Bits must be between 1 and 31");
        }
        let qmax = 1_u64 << bits;
        DimensionRange {
            continuous_minimum: min,
            continuous_maximum: max,
            quantized_maximum: qmax,
            scale: (qmax as f64) / (max - min),
            weight: weight
        }
    }

    pub fn to_quantized(&self, value: f64) -> Option<u32> {
        if value < self.continuous_minimum || value > self.continuous_maximum {
            None
        }
        else {
            let mut xformed = ((value - self.continuous_minimum) * self.scale).floor() as u64;
            if xformed >= self.quantized_maximum {
                xformed -= 1;
            }
            Some(xformed as u32)
        }
    }

    pub fn to_continuous(&self, value: u32) -> Option<f64> {
        if value as u64 >= self.quantized_maximum {
            None
        }
        else {
            Some(((value as f64) / self.scale) + self.continuous_minimum)
        }
    }

    /// Compute the weighted distance between two points 
    /// using weights derived from DimensionRanges.
    pub fn distance(point1: &[f64], point2: &[f64], ranges: &[DimensionRange]) -> f64 {
        let mut squared_diffs = 0_f64;
        let mut i = 0;
        for (a, b) in point1.iter().zip(point2) {
            let delta = ranges[i].weight * (a - b) * (a - b);
            squared_diffs += delta;
            i += 1;
        }
        squared_diffs.sqrt()
    }
}

pub struct Quantizer<F> where F: Fn(&[f64]) -> f64 {
    ranges: Vec<DimensionRange>,
    bits_per_dimension: usize,
    dimensions: usize,
    func: F
}

impl<F> Quantizer<F> where F: Fn(&[f64]) -> f64 {
    pub fn new(ranges: Vec<DimensionRange>, bits: usize, func: F) -> Self {
        let dimensions = ranges.len();
        Self {
            ranges: ranges,
            bits_per_dimension: bits,
            dimensions: dimensions,
            func: func
        }
    }
}

impl<F> HilbertModel for Quantizer<F> where F: Fn(&[f64]) -> f64 {
    /// Number of dimensions for points in the space.
    fn dimensions(&self) -> u16 {
        self.dimensions as u16
    }

    /// Number of bits of precision to use for each dimension.
    fn bits(&self) -> u8 {
        self.bits_per_dimension as u8
    }

    /// Map a point in the problem coordinate space to a Hilbert Index (via Hilbert Axes)
    fn to_hilbert_index(&self, point: &[f64]) -> num::BigUint {
        let mut axes = Vec::new();
        for (coord, range) in point.iter().zip(self.ranges.iter()) {
            let q = range.to_quantized(*coord).unwrap();
            axes.push(q);
        }
        // TODO: Store an interleaver in the Quantizer for resuse to improve performance.
        hilbert_index(&axes, self.bits() as usize, None)
    }

    /// Map a Hilbert Index to a point in problem coordinate space (via Hilbert Axes)
    fn to_point(&self, hilbert_index: &BigUint) -> Vec<f64> {
        let axes = hilbert_axes(hilbert_index, self.bits() as usize, self.dimensions() as usize);
        let mut point = Vec::new();
        for (q, range) in axes.iter().zip(self.ranges.iter()) {
            let coord = range.to_continuous(*q).unwrap();
            point.push(coord);
        }
        point
    }

    /// Evaluate the integrand function, returning NAN if the function is not defined for the given point.
    fn evaluate(&self, point: &[f64]) -> f64 {
        (self.func)(point)
    }

    /// Compute the distance between two points in the problem coordinate system, 
    /// which may be a weighted distsnce.
    fn distance(&self, point1: &[f64], point2: &[f64]) -> f64 {
        DimensionRange::distance(point1, point2, &self.ranges)
    }

    /// If the class implementing this trait cached any values, free them.
    fn clear(&self) {
        // No-op.
    }
}


#[cfg(test)]
mod tests {
    use super::DimensionRange;
    use super::Quantizer;
    use super::HilbertModel;

    #[test]
    fn hilbert_model_evaluate() {
        let model = create_test_model();
        let point = vec![3.0, 5.0];
        assert_eq!(model.evaluate(&point), 16.0);
    }

    #[test]
    fn hilbert_model_dimensions() {
        let model = create_test_model();
        assert_eq!(model.dimensions(), 2);
    }

    #[test]
    fn hilbert_model_bits() {
        let model = create_test_model();
        assert_eq!(model.bits(), 20);
    }

    #[test]
    fn hilbert_model_distance() {
        let model = create_test_model();
        let point1 = vec![3.0, 5.0];
        let point2 = vec![0.0, 1.0];
        assert_eq!(model.distance(&point1, &point2), 5.0);
    }

    #[test]
    fn dimension_range_to_quantized() {
        let range = DimensionRange::new(-10.0, 10.0, 1.0, 20);
        let actual_quantized = range.to_quantized(4.0).unwrap();
        let expected_quantized = 734003_u32; // ((14.0 / 20.0) * 2.0_f64.powi(20)).floor() as u32;
        assert_eq!(actual_quantized, expected_quantized);
    }

    #[test]
    fn dimension_range_to_continuous() {
        let range = DimensionRange::new(-10.0, 10.0, 1.0, 20);
        let quantized = 734003_u32;
        let expected_continuous = 4.0;
        let actual_continuous = range.to_continuous(quantized).unwrap();
        assert!((actual_continuous - expected_continuous).abs() < 0.00001);
    }

    /// Converting from Problem space to Hilbert index and back is lossy. 
    /// Verify that the deviation is small.
    #[test]
    fn hilbert_model_bidirectional() {
        let model = create_test_model();
        let point = vec![3.0, 5.0];
        let index = model.to_hilbert_index(&point);
        let actual_point = model.to_point(&index);
        let distance = model.distance(&point, &actual_point);
        
        assert!(distance.abs() < 0.00001);
    }

    fn create_test_model() -> Box<impl HilbertModel> {
        let f = |p: &[f64]| -> f64 { p[0]*p[1] + 1.0 };
        let ranges = vec![
            DimensionRange::new(-10.0, 10.0, 1.0, 20),
            DimensionRange::new(0.0, 10.0, 1.0, 20),
        ];
        Box::new(Quantizer::new(ranges, 20, f))
    }

}
