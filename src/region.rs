// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

use crate::common::{ib_to_bf_prec, PI};
use dashu_float::round::mode::{self, HalfEven};
use dashu_float::{Context, FBig};
use dashu_int::IBig;
use nalgebra::{Matrix2, Vector2};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::common::{fb_with_prec, get_prec_bits};
use std::fmt::{Debug, Display, Formatter};

/// Represents a closed interval [l, r] on the real line.
///
/// Used to constrain values in the GridSynth algorithm, particularly
/// in the ODGP (Orthogonal Diophantine Grid Problem) solver.
#[derive(Debug, Clone, PartialEq)]
pub struct Interval {
    pub l: FBig<HalfEven>,
    pub r: FBig<HalfEven>,
}

impl Interval {
    /// Creates a new interval [l, r].
    ///
    /// # Arguments
    ///
    /// * `l` - The left endpoint
    /// * `r` - The right endpoint
    pub fn new(l: FBig<HalfEven>, r: FBig<HalfEven>) -> Self {
        Self { l, r }
    }

    /// Computes the width of the interval.
    ///
    /// # Returns
    ///
    /// The difference r - l.
    pub fn width(&self) -> FBig<HalfEven> {
        &self.r - &self.l
    }

    /// Expands the interval by epsilon on both sides.
    ///
    /// # Arguments
    ///
    /// * `eps` - The amount to expand by
    ///
    /// # Returns
    ///
    /// A new interval [l - eps, r + eps].
    pub fn fatten(&self, eps: &FBig<HalfEven>) -> Self {
        Self {
            l: &self.l - eps,
            r: &self.r + eps,
        }
    }

    /// Checks if a value is within the interval.
    ///
    /// # Arguments
    ///
    /// * `x` - The value to check
    ///
    /// # Returns
    ///
    /// `true` if l ≤ x ≤ r, `false` otherwise.
    pub fn within(&self, x: &FBig<HalfEven>) -> bool {
        (self.l <= *x) && (*x <= self.r)
    }

    /// Scales the interval by a factor.
    ///
    /// # Arguments
    ///
    /// * `factor` - The scaling factor
    ///
    /// # Returns
    ///
    /// A new interval [l*factor, r*factor] if factor ≥ 0,
    /// or [r*factor, l*factor] if factor < 0.
    pub fn scale(&self, factor: &FBig<HalfEven>) -> Self {
        let zero = ib_to_bf_prec(IBig::ZERO);
        if *factor >= zero {
            Self {
                l: &self.l * factor,
                r: &self.r * factor,
            }
        } else {
            Self {
                l: &self.r * factor,
                r: &self.l * factor,
            }
        }
    }

    /// Subtracts a value from the interval.
    ///
    /// # Arguments
    ///
    /// * `x` - The value to subtract
    ///
    /// # Returns
    ///
    /// A new interval [l - x, r - x].
    pub fn sub_ref(&self, x: &FBig<HalfEven>) -> Self {
        Interval {
            l: &self.l - x,
            r: &self.r - x,
        }
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{}, {}]", self.l, self.r)
    }
}

// Interval + Interval
impl Add for Interval {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            l: &self.l + &rhs.l,
            r: &self.r + &rhs.r,
        }
    }
}

// Interval + f64
impl Add<f64> for Interval {
    type Output = Self;
    fn add(self, rhs: f64) -> Self {
        let rhs_float = FBig::<HalfEven>::try_from(rhs).unwrap();
        Self {
            l: &self.l + &rhs_float,
            r: &self.r + &rhs_float,
        }
    }
}

// Interval + FBig
impl Add<FBig<HalfEven>> for Interval {
    type Output = Self;
    fn add(self, rhs: FBig<HalfEven>) -> Self {
        Self {
            l: &self.l + &rhs,
            r: &self.r + &rhs,
        }
    }
}

// f64 + Interval
impl Add<Interval> for f64 {
    type Output = Interval;
    fn add(self, rhs: Interval) -> Interval {
        rhs + self
    }
}

// FBig + Interval
impl Add<Interval> for FBig<HalfEven> {
    type Output = Interval;
    fn add(self, rhs: Interval) -> Interval {
        rhs + self
    }
}

// Interval - Interval（= Interval + (-Interval)）
impl Sub for Interval {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

// f64 - Interval
impl Sub<Interval> for f64 {
    type Output = Interval;
    fn sub(self, rhs: Interval) -> Interval {
        -rhs + self
    }
}

// FBig - Interval
impl Sub<Interval> for FBig<HalfEven> {
    type Output = Interval;
    fn sub(self, rhs: Interval) -> Interval {
        -rhs + self
    }
}

// Interval - f64
impl Sub<f64> for Interval {
    type Output = Interval;
    fn sub(self, rhs: f64) -> Interval {
        let rhs_float = FBig::<HalfEven>::try_from(rhs).unwrap();
        self + (-Interval::new(rhs_float.clone(), rhs_float))
    }
}

// Interval - IBig
impl Sub<IBig> for Interval {
    type Output = Interval;
    fn sub(self, rhs: IBig) -> Interval {
        self + (-Interval::new(ib_to_bf_prec(rhs.clone()), ib_to_bf_prec(rhs)))
    }
}

// Interval - FBig
impl Sub<FBig<HalfEven>> for Interval {
    type Output = Interval;
    fn sub(self, rhs: FBig<HalfEven>) -> Interval {
        self + (-Interval::new(rhs.clone(), rhs))
    }
}

// -Interval
impl Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            l: -self.r,
            r: -self.l,
        }
    }
}

// Interval * f64
impl Mul<f64> for Interval {
    type Output = Interval;
    fn mul(self, rhs: f64) -> Interval {
        if rhs >= 0.0 {
            Interval {
                l: self.l * FBig::<HalfEven>::try_from(rhs).unwrap(),
                r: self.r * FBig::<HalfEven>::try_from(rhs).unwrap(),
            }
        } else {
            Interval {
                l: self.r * FBig::<HalfEven>::try_from(rhs).unwrap(),
                r: self.l * FBig::<HalfEven>::try_from(rhs).unwrap(),
            }
        }
    }
}

// Interval * FBig
impl Mul<FBig<HalfEven>> for Interval {
    type Output = Interval;
    fn mul(self, rhs: FBig<HalfEven>) -> Interval {
        if rhs >= ib_to_bf_prec(IBig::ZERO) {
            Interval {
                l: self.l * rhs.clone(),
                r: self.r * rhs,
            }
        } else {
            Interval {
                l: self.r * rhs.clone(),
                r: self.l * rhs,
            }
        }
    }
}

// f64 * Interval
impl Mul<Interval> for f64 {
    type Output = Interval;
    fn mul(self, rhs: Interval) -> Interval {
        rhs * self
    }
}

// FBig * Interval
impl Mul<Interval> for FBig<HalfEven> {
    type Output = Interval;
    fn mul(self, rhs: Interval) -> Interval {
        rhs * self
    }
}

// Interval / FBig
impl Div<FBig<HalfEven>> for Interval {
    type Output = Interval;
    fn div(self, rhs: FBig<HalfEven>) -> Interval {
        if rhs > ib_to_bf_prec(IBig::ZERO) {
            Interval {
                l: self.l / rhs.clone(),
                r: self.r / rhs,
            }
        } else {
            Interval {
                l: self.r / rhs.clone(),
                r: self.l / rhs,
            }
        }
    }
}

/// Represents a rectangle in the 2D plane as a product of two intervals.
///
/// Used to represent bounding boxes and constraint regions in the GridSynth algorithm.
#[derive(Debug, Clone, PartialEq)]
pub struct Rectangle {
    pub x: Interval,
    pub y: Interval,
}

impl Rectangle {
    /// Creates a new rectangle from interval endpoints.
    ///
    /// # Arguments
    ///
    /// * `x_l` - Left x-coordinate
    /// * `x_r` - Right x-coordinate
    /// * `y_l` - Bottom y-coordinate
    /// * `y_r` - Top y-coordinate
    pub fn new(
        x_l: FBig<HalfEven>,
        x_r: FBig<HalfEven>,
        y_l: FBig<HalfEven>,
        y_r: FBig<HalfEven>,
    ) -> Self {
        Self {
            x: Interval::new(x_l, x_r),
            y: Interval::new(y_l, y_r),
        }
    }

    /// Computes the area of the rectangle.
    ///
    /// # Returns
    ///
    /// The product of the widths of the x and y intervals.
    /// Computes the area of the ellipse.
    ///
    /// # Returns
    ///
    /// π / √(det(d)).
    pub fn area(&self) -> FBig<HalfEven> {
        self.x.width() * self.y.width()
    }

    fn scale(&self, factor: FBig<HalfEven>) -> Self {
        if factor >= ib_to_bf_prec(IBig::ZERO) {
            Self {
                x: self.x.scale(&factor),
                y: self.y.scale(&factor),
            }
        } else {
            Self {
                x: self.y.scale(&factor),
                y: self.x.scale(&factor),
            }
        }
    }
}

impl fmt::Display for Rectangle {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}×{}", self.x, self.y)
    }
}

// Rectangle * FBig
impl Mul<FBig<HalfEven>> for Rectangle {
    type Output = Self;

    fn mul(self, rhs: FBig<HalfEven>) -> Self::Output {
        self.scale(rhs)
    }
}

// FBig * Rectangle
impl Mul<Rectangle> for FBig<HalfEven> {
    type Output = Rectangle;

    fn mul(self, rhs: Rectangle) -> Rectangle {
        rhs.scale(self)
    }
}

/// Represents an ellipse in the 2D plane.
///
/// An ellipse is defined by a positive definite matrix d and a center point p.
/// Points (x, y) are inside the ellipse if (v - p)ᵀ d (v - p) ≤ 1.
#[derive(Debug, Clone)]
pub struct Ellipse {
    pub d: Matrix2<FBig<HalfEven>>,
    pub p: Vector2<FBig<HalfEven>>,
}

impl Ellipse {
    /// Creates a new ellipse from a matrix and center point.
    ///
    /// # Arguments
    ///
    /// * `d` - The positive definite matrix defining the ellipse shape
    /// * `p` - The center point of the ellipse
    pub fn new(d: Matrix2<FBig<HalfEven>>, p: Vector2<FBig<HalfEven>>) -> Self {
        Self { d, p }
    }

    /// Creates a new ellipse from individual matrix and point components.
    ///
    /// # Arguments
    ///
    /// * `d00`, `d01`, `d10`, `d11` - Matrix elements
    /// * `px`, `py` - Center point coordinates
    pub fn from(
        d00: FBig<HalfEven>,
        d01: FBig<HalfEven>,
        d10: FBig<HalfEven>,
        d11: FBig<HalfEven>,
        px: FBig<HalfEven>,
        py: FBig<HalfEven>,
    ) -> Self {
        Self {
            d: Matrix2::new(d00, d01, d10, d11),
            p: Vector2::new(px, py),
        }
    }

    pub fn px(&self) -> &FBig<HalfEven> {
        &self.p[0]
    }

    pub fn py(&self) -> &FBig<HalfEven> {
        &self.p[1]
    }

    pub fn a(&self) -> &FBig<HalfEven> {
        &self.d[(0, 0)]
    }

    pub fn b(&self) -> &FBig<HalfEven> {
        &self.d[(0, 1)]
    }

    pub fn d(&self) -> &FBig<HalfEven> {
        &self.d[(1, 1)]
    }

    /// Computes the skew of the ellipse.
    ///
    /// # Returns
    ///
    /// The square of the off-diagonal element b².
    pub fn skew(&self) -> FBig<HalfEven> {
        self.b().powi(IBig::from(2))
    }

    /// Computes the bias (aspect ratio) of the ellipse.
    ///
    /// # Returns
    ///
    /// The ratio d/a of the diagonal elements.
    pub fn bias(&self) -> FBig<HalfEven> {
        self.d() / self.a()
    }
    /// Checks if a point is inside the ellipse.
    ///
    /// # Arguments
    ///
    /// * `v` - The point to check
    ///
    /// # Returns
    ///
    /// `true` if (v - p)ᵀ d (v - p) ≤ 1, `false` otherwise.
    pub fn inside(&self, v: &Vector2<FBig<HalfEven>>) -> bool {
        let x = &v[0] - &self.p[0];
        let y = &v[1] - &self.p[1];
        let a = self.a();
        let b = self.b();
        let d = self.d();
        let x_sq = &x * &x;
        let term1 = a * &x_sq;

        let x_y = &x * &y;
        let two_b = 2 * b;
        let term2 = &two_b * &x_y;

        let y_sq = &y * &y;
        let term3 = d * &y_sq;

        let sum12 = &term1 + &term2;
        let value = &sum12 + &term3;
        value <= fb_with_prec(FBig::<HalfEven>::from(1))
    }

    /// Computes the axis-aligned bounding box of the ellipse.
    ///
    /// # Returns
    ///
    /// A `Rectangle` that tightly bounds the ellipse.
    pub fn bbox(&self) -> Rectangle {
        let ctx: Context<mode::HalfEven> = Context::<mode::HalfEven>::new(get_prec_bits());
        let sqrt_det = self.sqrt_det();
        let w = ctx.sqrt(self.d().repr()).value() / &sqrt_det;
        let h = ctx.sqrt(self.a().repr()).value() / &sqrt_det;
        let px_minus_w = self.px() - &w;
        let px_plus_w = self.px() + &w;
        let py_minus_h = self.py() - &h;
        let py_plus_h = self.py() + &h;
        Rectangle {
            x: Interval::new(px_minus_w, px_plus_w),
            y: Interval::new(py_minus_h, py_plus_h),
        }
    }

    /// Computes the square root of the determinant of the matrix d.
    ///
    /// # Returns
    ///
    /// √(det(d)) = √(ad - b²).
    pub fn sqrt_det(&self) -> FBig<HalfEven> {
        let ctx: Context<mode::HalfEven> = Context::<mode::HalfEven>::new(get_prec_bits());
        let det = self.d() * self.a() - self.b().powi(IBig::from(2));
        ctx.sqrt(det.repr()).value()
    }

    pub fn area(&self) -> FBig<HalfEven> {
        (*PI).clone() / self.sqrt_det()
    }

    /// Normalizes the ellipse to have unit area.
    ///
    /// # Returns
    ///
    /// A new ellipse with the same shape but scaled to have area π.
    pub fn normalize(&self) -> Self {
        let ctx: Context<mode::HalfEven> = Context::<mode::HalfEven>::new(get_prec_bits());
        let factor = self.sqrt_det();
        let factor_sqrt = ctx.sqrt(factor.repr()).value();
        Ellipse::new(self.d.clone() / factor, self.p.clone() * factor_sqrt)
    }
}

// Ellipse * FBig
impl Mul<FBig<HalfEven>> for Ellipse {
    type Output = Ellipse;
    fn mul(self, rhs: FBig<HalfEven>) -> Ellipse {
        let inv_rhs: FBig<HalfEven> = 1 / &rhs;
        let inv_rhs_sq = inv_rhs.powi(IBig::from(2));
        Ellipse::new(self.d * inv_rhs_sq, self.p * rhs)
    }
}

// FBig * Ellipse
impl Mul<Ellipse> for FBig<HalfEven> {
    type Output = Ellipse;
    fn mul(self, rhs: Ellipse) -> Ellipse {
        rhs * self
    }
}

// Ellipse / FBig
impl Div<FBig<HalfEven>> for Ellipse {
    type Output = Ellipse;
    fn div(self, rhs: FBig<HalfEven>) -> Ellipse {
        let rhs_sq = rhs.clone().powi(IBig::from(2));
        Ellipse::new(self.d * rhs_sq, self.p / rhs)
    }
}

impl Display for Ellipse {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Ellipse {{\n  d: [[{}, {}],\n      [{}, {}]],\n  p: ({}, {})\n}}",
            self.d[(0, 0)].to_decimal().value(),
            self.d[(0, 1)].to_decimal().value(),
            self.d[(1, 0)].to_decimal().value(),
            self.d[(1, 1)].to_decimal().value(),
            self.p[0].to_decimal().value(),
            self.p[1].to_decimal().value(),
        )
    }
}
