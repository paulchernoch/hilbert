use criterion::{black_box, criterion_group, criterion_main, Criterion};
use num::BigUint;
use rand::thread_rng;
use rand::seq::SliceRandom;
use hilbert::transform::fast_hilbert;
use hilbert::interleaver::Interleaver;
use hilbert::point::Point;

fn point_square_distance(p1 : &Point, p2 : &Point) -> u64 {
    p1.square_distance_loop_unrolling(p2)
}

fn hilbert_round_trip(index : u32, dimensions : usize, bits : usize, interleaver_option : Option<&Interleaver>) {
    let hilbert_index_in : BigUint = index.into();
    let coordinates = fast_hilbert::hilbert_axes(&hilbert_index_in, bits, dimensions);
    let _hilbert_index_out = fast_hilbert::hilbert_index(&coordinates, bits, interleaver_option);
}

fn test_points(n : usize, dimensions : u32) -> Vec<Point> {
    let mut rng = thread_rng();
    let points : Vec<Point> = (0..n).map(|i| {
        let mut point_data : Vec<u32> = (0..dimensions).collect();
        point_data.shuffle(&mut rng);
        Point::new(i, &point_data)
    } ).collect();
    points
} 

fn criterion_benchmark(c: &mut Criterion) {
    let points = test_points(2, 100);
    let (p1, p2) = (&points[0], &points[1]);

    c.bench_function("point square distance, 100 dimensions", |b| b.iter(|| point_square_distance(p1, p2)));

    // Switch between true and false between benchmark runs to see the impact.
    let use_interleaver = false;

    let interleaver_50_10 = Interleaver::new(50, 10);
    let interleaver_100_10 = Interleaver::new(100, 10);
    let interleaver_200_10 = Interleaver::new(200, 10);

    {
        let interleaver_50_10_opt = if use_interleaver { Some(&interleaver_50_10) } else { None };
        let interleaver_100_10_opt = if use_interleaver { Some(&interleaver_100_10) } else { None };
        let interleaver_200_10_opt = if use_interleaver { Some(&interleaver_200_10) } else { None };

        c.bench_function("hilbert 50 D 10 B", |b| b.iter(|| hilbert_round_trip(black_box(12345), 50, 10, interleaver_50_10_opt)));
        c.bench_function("hilbert 100 D 10 B", |b| b.iter(|| hilbert_round_trip(black_box(12345), 100, 10, interleaver_100_10_opt)));
        c.bench_function("hilbert 200 D 10 B", |b| b.iter(|| hilbert_round_trip(black_box(12345), 200, 10, interleaver_200_10_opt)));
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
