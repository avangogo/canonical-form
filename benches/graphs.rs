use canonical_form::Canonize;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

#[derive(Ord, PartialOrd, PartialEq, Eq, Clone, Debug)]
struct Graph {
    adj: Vec<Vec<usize>>,
}

impl Graph {
    fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut adj = vec![Vec::new(); n];
        for &(u, v) in edges {
            adj[u].push(v);
            adj[v].push(u);
        }
        Graph { adj }
    }
}

impl Canonize for Graph {
    fn size(&self) -> usize {
        self.adj.len()
    }
    fn apply_morphism(&self, perm: &[usize]) -> Self {
        let mut adj = vec![Vec::new(); self.size()];
        for (i, nbrs) in self.adj.iter().enumerate() {
            adj[perm[i]] = nbrs.iter().map(|&u| perm[u]).collect();
            adj[perm[i]].sort();
        }
        Graph { adj }
    }
    fn invariant_neighborhood(&self, u: usize) -> Vec<Vec<usize>> {
        vec![self.adj[u].clone()]
    }
}

fn petersen() -> Graph {
    Graph::new(
        10,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 0),
            (5, 7),
            (7, 9),
            (9, 6),
            (6, 8),
            (8, 5),
            (0, 5),
            (1, 6),
            (2, 7),
            (3, 8),
            (4, 9),
        ],
    )
}

fn empty(n: usize) -> Graph {
    Graph::new(n, &[])
}

fn clique(n: usize) -> Graph {
    let mut edges = Vec::with_capacity((n * (n - 1)) / 2);
    for i in 0..n {
        for j in 0..i {
            edges.push((i, j))
        }
    }
    Graph::new(n, &edges)
}

fn cycle(n: usize) -> Graph {
    let mut edges = Vec::with_capacity(n);
    for i in 1..n {
        edges.push((i - 1, i))
    }
    edges.push((n - 1, 0));
    Graph::new(n, &edges)
}

extern crate rand;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn gnp<R: Rng>(n: usize, p: f64, rng: &mut R) -> Graph {
    let mut edges = Vec::new();
    for i in 0..n {
        for j in 0..i {
            if rng.gen_bool(p) {
                edges.push((i, j))
            }
        }
    }
    Graph::new(n, &edges)
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let petersen = petersen();
    c.bench_function("petersen", |b| b.iter(|| black_box(&petersen).canonical()));
    let empty = empty(20);
    c.bench_function("empty 20", |b| b.iter(|| black_box(&empty).canonical()));
    let clique_ = clique(20);
    c.bench_function("clique 20", |b| b.iter(|| black_box(&clique_).canonical()));
    let cycle = cycle(50);
    c.bench_function("cycle 50", |b| b.iter(|| black_box(&cycle).canonical()));
    let mut rng: StdRng = SeedableRng::from_seed([42; 32]);
    let gnp_ = gnp(200, 0.1, &mut rng);
    c.bench_function("gnp 200", |b| b.iter(|| black_box(&gnp_).canonical()));
    let graphs7: Vec<_> = (0..100).map(|_| gnp(10, 0.5, &mut rng)).collect();
    c.bench_function("graphs 10", |b| {
        b.iter(|| {
            for g in &graphs7 {
                black_box(g).canonical();
            }
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
