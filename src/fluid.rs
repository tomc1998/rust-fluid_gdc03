#![allow(dead_code)]

/// Macro for indexing into a 1D array using 2D coordinates.
macro_rules! IX {
  ( $x: expr, $y: expr,  $w: expr ) => {
    { $x as usize + $y as usize * $w as usize }
  }
}

/// A source of fluid.
pub struct Source {
  /// The index of the source in the grid array
  pub ix: usize,
  /// Value between 0 and 1, density of the fluid at source.
  pub density: f32,
}

impl Source {
  /// Create a new source, given an index into a float array for its position.
  pub fn new(ix: usize, density: f32) -> Source {
    Source { ix: ix, density: density }
  }
}

/// Add borders to the density grid. We can do this by just setting the density
/// of the outer cells to be equal to the density of the cells inset by 1. This
/// means the difference between the two cells is 0, so there will never be any
/// flow outwards.
/// # Params
/// * `b` - The type of border. 1 for vertical vel walls, 2 for hori vel walls, 0 for dens.
fn set_borders(grid: &mut [f32], grid_w: u32, b: u8) {
  let grid_h = grid.len()/grid_w as usize;
  for ii in 0..grid_w {
    let ix_top = IX!(ii, 0, grid_w);
    let ix_top_inset = IX!(ii, 1, grid_w);
    let ix_bot = IX!(ii, grid_h-1, grid_w);
    let ix_bot_inset = IX!(ii, grid_h-2, grid_w);
    grid[ix_top] = {if b == 2 {-grid[ix_top_inset]} else {grid[ix_top_inset]}};
    grid[ix_bot] = {if b == 2 {-grid[ix_bot_inset]} else {grid[ix_bot_inset]}};
  }
  for ii in 0..grid_h {
    let ix_left = IX!(0, ii, grid_w);
    let ix_left_inset = IX!(1, ii, grid_w);
    let ix_right = IX!(grid_w-1, ii, grid_w);
    let ix_right_inset = IX!(grid_w-2, ii, grid_w);
    grid[ix_left] = {if b == 1 {-grid[ix_left_inset]} else {grid[ix_left_inset]}};
    grid[ix_right] = {if b == 1 {-grid[ix_right_inset]} else {grid[ix_right_inset]}};
  }
  grid[IX!(0, 0, grid_w)]             = 0.5*(grid[IX!(1, 0, grid_w)]+grid[IX!(0, 1, grid_w)]);
  grid[IX!(0, grid_h-1, grid_w)]      = 0.5*(grid[IX!(1, grid_h-1, grid_w)]+grid[IX!(0, grid_h-2, grid_w)]);
  grid[IX!(grid_w-1, 0, grid_w)]      = 0.5*(grid[IX!(grid_w-2, 0, grid_w)]+grid[IX!(grid_w-1, 1, grid_w)]);
  grid[IX!(grid_w-1, grid_h-1, grid_w)] = 0.5*(grid[IX!(grid_w-2, grid_h-1, grid_w)]+grid[IX!(grid_w-1, grid_h-2, grid_w)]);
}

/// Process diffusion
/// # Params
/// * `b` - Border type, see set_borders
fn diffuse(grid: &mut [f32], prev_grid: &mut [f32], grid_w: u32, dt: f32, diff: f32, borders: bool, b: u8) {
  let a = dt*diff*grid.len() as f32;
  let grid_h = grid.len()/grid_w as usize;
  let iterations = 1;
  if diff == 0.0 { return; }
  // Iteratively diffuse, solve the linear system to find the next state of the grid
  for _ in 0..iterations {
    for ii in 1..grid_w-1 {
      for jj in 1..grid_h-1 {
        let ix = IX!(ii, jj, grid_w);
        let ix_up = IX!(ii, jj-1, grid_w); // 1 row up
        let ix_down = IX!(ii, jj+1, grid_w); // 1 row down
        grid[ix] = (prev_grid[ix] + a*(grid[ix-1] + grid[ix+1] 
                                       + grid[ix_up] + grid[ix_down]))/(1.0+4.0*a);
      }
    }
    if borders { set_borders(grid, grid_w, b) }
  }
}

/// Process density movement via velocity
/// # Params
/// * `b` - Border type, see set_borders
fn advect(grid: &mut [f32], prev_grid: &[f32], vx_grid: &[f32], vy_grid: &[f32], grid_w: u32, dt: f32, borders: bool, b: u8) {
  let grid_h = grid.len()/grid_w as usize;
  let dt0 = dt*grid_w as f32;

  for ii in 1..grid_w-1 {
    for jj in 1..grid_h-1 {
      let ix = IX!(ii, jj, grid_w);
      let mut x = ii as f32 - dt0 * vx_grid[ix];
      let mut y = jj as f32 - dt0 * vy_grid[ix];

      if x < 0.5 { x = 0.5 }
      if x > grid_w as f32 - 1.5 { x = grid_w as f32 - 1.5 }
      let ii0 = x as u32;
      let ii1 = ii0 + 1;
      if y < 0.5 { y = 0.5 }
      if y > grid_h as f32 - 1.5 { y = grid_h as f32 - 1.5 }
      let jj0 = y as u32;
      let jj1 = jj0 + 1;

      let s1 = x - ii0 as f32;
      let s0 = 1.0 - s1;
      let t1 = y - jj0 as f32;
      let t0 = 1.0 - t1;

      let ix00 = IX!(ii0, jj0, grid_w);
      let ix01 = IX!(ii0, jj1, grid_w);
      let ix11 = IX!(ii1, jj1, grid_w);
      let ix10 = IX!(ii1, jj0, grid_w);
      grid[ix] = s0 * (t0*prev_grid[ix00]+t1*prev_grid[ix01])
        + s1 * (t0*prev_grid[ix10] + t1*prev_grid[ix11]);
    }
  }
  if borders { set_borders(grid, grid_w, b); }
}

fn project(vx_grid: &mut [f32], vy_grid: &mut [f32], p: &mut [f32], div: &mut [f32], grid_w: u32, borders: bool) {
  let h = 1.0/grid_w as f32;
  let grid_h = vx_grid.len()/grid_w as usize;
  for ii in 1..grid_w-1 {
    for jj in 1..grid_h-1 {
      let ix = IX!(ii, jj, grid_w);
      let ix_up = IX!(ii, jj-1, grid_w); // 1 row up
      let ix_down = IX!(ii, jj+1, grid_w); // 1 row down
      div[ix] = -0.5*h*(vx_grid[ix+1]-vx_grid[ix-1] +
                        vy_grid[ix_down]-vy_grid[ix_up]);
      p[ix] = 0.0;
    }
  }
  if borders { set_borders(div, grid_w, 0); set_borders(p, grid_w, 0) }

  let iterations = 1;
  for _ in 0..iterations {
    for ii in 1..grid_w-1 {
      for jj in 1..grid_h-1 {
        let ix = IX!(ii, jj, grid_w);
        let ix_up = IX!(ii, jj-1, grid_w); // 1 row up
        let ix_down = IX!(ii, jj+1, grid_w); // 1 row down
        p[ix] = (div[ix] + p[ix-1] + p[ix+1] + p[ix_up] + p[ix_down])/4.0;
      }
    }
    if borders { set_borders(p, grid_w, 0) }
  }

  for ii in 1..grid_w-1 {
    for jj in 1..grid_h-1 {
      let ix = IX!(ii, jj, grid_w);
      let ix_up = IX!(ii, jj-1, grid_w); // 1 row up
      let ix_down = IX!(ii, jj+1, grid_w); // 1 row down
      vx_grid[ix] -= 0.5*(p[ix+1]-p[ix-1])/h;
      vy_grid[ix] -= 0.5*(p[ix_down]-p[ix_up])/h;
    }
  }
  if borders { set_borders(vx_grid, grid_w, 1); set_borders(vy_grid, grid_w, 2) }
}

/// Step density
fn step_dens(dens_grid: &mut [f32], vx_grid: &[f32], vy_grid: &[f32], grid_w: u32, dt: f32, diff: f32, borders: bool) {
  // Make a copy of the dens_grid
  let mut prev_dens_grid = dens_grid.to_vec();
  let prev_dens_grid = &mut prev_dens_grid[..];

  // Swap binding in preparation for the next swap
  //let (prev_dens_grid, dens_grid) = (dens_grid, prev_dens_grid);

  // Process diffusion
  diffuse(dens_grid, prev_dens_grid, grid_w, dt, diff, borders, 0);

  // Swap bindings, b/c we just updates dens_grid and advect() needs to use
  // that for the previous grid
  //let (prev_dens_grid, dens_grid) = (dens_grid, prev_dens_grid);

  // Process velocity of particles
  advect(dens_grid, prev_dens_grid, vx_grid, vy_grid, grid_w, dt, borders, 0);
}

/// Step velocity
fn step_vel(vx_grid: &mut [f32], vy_grid: &mut [f32], grid_w: u32, dt: f32, diff: f32, borders: bool) {
  // Copy velocity grids, so we have a copy of it before processing each step
  let mut prev_vx_grid = vx_grid.to_vec();
  let mut prev_vy_grid = vy_grid.to_vec();
  let prev_vx_grid = &mut prev_vx_grid[..];
  let prev_vy_grid = &mut prev_vy_grid[..];

  let (prev_vx_grid, vx_grid) = (vx_grid, prev_vx_grid);
  let (prev_vy_grid, vy_grid) = (vy_grid, prev_vy_grid);

  // Diffuse just like with density but with velocity instead
  //diffuse(vx_grid, prev_vx_grid, grid_w, dt, diff, borders, 1);
  //diffuse(vy_grid, prev_vy_grid, grid_w, dt, diff, borders, 2);
  project(vx_grid, vy_grid, prev_vx_grid, prev_vy_grid, grid_w, borders);

  let (prev_vx_grid, vx_grid) = (vx_grid, prev_vx_grid);
  let (prev_vy_grid, vy_grid) = (vy_grid, prev_vy_grid);

  // Advect just like with density
  advect(vx_grid, prev_vx_grid, prev_vx_grid, prev_vy_grid, grid_w, dt, borders, 1);
  advect(vy_grid, prev_vy_grid, prev_vx_grid, prev_vy_grid, grid_w, dt, borders, 2);
  project(vx_grid, vy_grid, prev_vx_grid, prev_vy_grid, grid_w, borders);
}

/// Steps a grid of floats containing fluid data according to dt.
pub fn step_fluid(dens_grid: &mut [f32], vx_grid: &mut [f32], vy_grid: &mut [f32], grid_w: u32, dt: f32, diff: f32, borders: bool) {
  // Step density, alter density grid
  step_dens(dens_grid, vx_grid, vy_grid, grid_w, dt, diff, borders);
  step_vel(vx_grid, vy_grid, grid_w, dt, diff, borders);
}
