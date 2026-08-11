// Copyright (c) 2024-2025 Shun Yamamoto and Nobuyuki Yoshioka, and IBM
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

pub mod fallback;
pub mod mixed_diagonal;
pub mod mixed_fallback;
pub mod mixing;

pub use crate::accuracy::AchievedDiamondError;
pub use mixing::{
    achieved_diagonal_diamond_error, diagonal_diamond_distance, diamond_to_spec_epsilon,
    mixture_weight, pauli_diamond_distance, MixtureWeight, WFrame,
};

pub use mixed_diagonal::{
    synth_mixed_diagonal, MixedDiagonalBranch, MixedDiagonalRegion, MixedDiagonalResult,
};

pub use fallback::{exact_q, synth_fallback, FallbackResult, SectorRegion};

pub use mixed_fallback::{synth_mixed_fallback, MixedFallbackResult, MixedFallbackSide};
