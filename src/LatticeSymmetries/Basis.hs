module LatticeSymmetries.Basis where

import Foreign.C.Types
import Foreign.Marshal.Utils (fromBool)

data SpinBasis = SpinBasis
  { sbNumberSites :: !Int,
    sbMagnetization :: !(Maybe Int)
  }

data SpinfulOccupation
  = SpinfulNoOccupation
  | SpinfulTotalParticles !Int
  | SpinfulPerSector !Int !Int

data SpinfulFermionicBasis = SpinfulFermionicBasis
  { sfbNumberSites :: !Int,
    sfbOccupation :: !SpinfulOccupation
  }

data SpinlessFermionicBasis = SpinlessFermionicBasis
  { fbNumberSites :: !Int,
    fbOccupation :: !(Maybe Int)
  }

newtype Cparticle_type = Cparticle_type CInt

c_LS_HS_SPIN :: Cparticle_type
c_LS_HS_SPIN = Cparticle_type 0

c_LS_HS_FERMION :: Cparticle_type
c_LS_HS_FERMION = Cparticle_type 1

data Cbasis = Cbasis
  { cbasis_number_sites :: {-# UNPACK #-} !CInt,
    cbasis_number_particles :: {-# UNPACK #-} !CInt,
    cbasis_number_up :: {-# UNPACK #-} !CInt,
    cbasis_particle_type :: {-# UNPACK #-} !Cparticle_type,
    cbasis_state_index_is_identity :: {-# UNPACK #-} !CBool
  }

class IsBasis a where
  toCbasis :: a -> Cbasis

instance IsBasis SpinBasis where
  toCbasis (SpinBasis n (Just m)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral n,
        -- we have: cbasis_number_up + cbasis_number_down = n
        --          cbasis_number_up - cbasis_number_down = m
        cbasis_number_up = fromIntegral $ (n + m) `div` 2,
        cbasis_particle_type = c_LS_HS_SPIN,
        cbasis_state_index_is_identity = fromBool False
      }
  toCbasis (SpinBasis n Nothing) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral n,
        cbasis_number_up = -1,
        cbasis_particle_type = c_LS_HS_SPIN,
        cbasis_state_index_is_identity = fromBool True
      }

instance IsBasis SpinfulFermionicBasis where
  toCbasis (SpinfulFermionicBasis n SpinfulNoOccupation) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = -1,
        cbasis_number_up = -1,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool True
      }
  toCbasis (SpinfulFermionicBasis n (SpinfulTotalParticles p)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral p,
        cbasis_number_up = -1,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool False
      }
  toCbasis (SpinfulFermionicBasis n (SpinfulPerSector u d)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral (u + d),
        cbasis_number_up = fromIntegral u,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool False
      }

instance IsBasis SpinlessFermionicBasis where
  toCbasis (SpinlessFermionicBasis n (Just p)) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = fromIntegral p,
        cbasis_number_up = 0,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool False
      }
  toCbasis (SpinlessFermionicBasis n Nothing) =
    Cbasis
      { cbasis_number_sites = fromIntegral n,
        cbasis_number_particles = -1,
        cbasis_number_up = 0,
        cbasis_particle_type = c_LS_HS_FERMION,
        cbasis_state_index_is_identity = fromBool True
      }
