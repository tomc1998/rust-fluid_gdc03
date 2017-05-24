#[macro_use]
extern crate glium;

mod fluid;

use glium::backend::glutin_backend::GlutinFacade;

#[derive(Copy, Clone)]
struct Vertex {
    pos: [f32; 2],
    uv: [f32; 2],
}
implement_vertex!(Vertex, pos, uv);

static VERT_SRC : &'static str = r#"
    #version 140

    in vec2 pos;
    in vec2 uv;

    out vec2 v_uv;

    void main() {
      v_uv = uv;
      gl_Position = vec4(pos, 0.0, 1.0);
    }
"#;

static FRAG_SRC : &'static str = r#"
    #version 140

    uniform highp sampler2D tex;
    uniform highp sampler2D tex1;
    uniform highp sampler2D tex2;

    in vec2 v_uv;

    out vec4 color;

    void main() {
      color = vec4(texture(tex, v_uv).x, texture(tex1, v_uv).x, texture(tex2, v_uv).x, 1.0);
    }
"#;


fn setup_display() -> GlutinFacade {
  use glium::DisplayBuild;
  glium::glutin::WindowBuilder::new().build_glium().unwrap()
}

fn setup_shader(display: &GlutinFacade) -> glium::Program {
  glium::Program::from_source(display, VERT_SRC, FRAG_SRC, None).unwrap()
}

fn main() {
  let display = setup_display();
  let shader = setup_shader(&display);
  let (display_w, display_h) = display.get_window().unwrap().get_inner_size().unwrap();

  let vbo_data = vec![
    Vertex{ pos: [-1.0, -1.0], uv: [0.0, 0.0] },
    Vertex{ pos: [ 1.0, -1.0], uv: [1.0, 0.0] },
    Vertex{ pos: [ 1.0,  1.0], uv: [1.0, 1.0] },
    Vertex{ pos: [-1.0, -1.0], uv: [0.0, 0.0] },
    Vertex{ pos: [-1.0,  1.0], uv: [0.0, 1.0] },
    Vertex{ pos: [ 1.0,  1.0], uv: [1.0, 1.0] }, ];


  let vbo = glium::VertexBuffer::new(&display, &vbo_data).unwrap();
  let indices = glium::index::NoIndices(glium::index::PrimitiveType::TrianglesList);


  // Set up fluids

  const TEX_SIZE: usize = 256;
  // Velocity grids
  let mut vx_grid = [0.0; TEX_SIZE*TEX_SIZE].to_vec();
  let mut vy_grid = [0.0; TEX_SIZE*TEX_SIZE].to_vec();
  let mut tex_data = [0.0; TEX_SIZE*TEX_SIZE].to_vec();
  // Create texture data buffer for fluid

  let (mut prev_mx, mut prev_my) = (0, 0);
  // 4 tuple, ABXY, AB for position and XY for velocity. Every frame the
  //   velocity cell this correcponds to gets set to XY.
  let mut curr_vel = (0, 0, 0.0, 0.0);
  loop {
    // listing the events produced by the window and waiting to be received
    for ev in display.poll_events() {
      match ev {
        glium::glutin::Event::Closed => { return }   // the window has been closed by the user
        glium::glutin::Event::MouseMoved(x, y) => {
          let x = (x as f32 / display_w as f32)*TEX_SIZE as f32;
          let y = ((display_h as f32 - y as f32) / display_h as f32)*TEX_SIZE as f32;
          let (x, y) = (x as usize, y as usize);
          let ix = x + y*TEX_SIZE;
          if ix >= tex_data.len() { continue; }
          tex_data[x + y*TEX_SIZE] = 1.0;
          if prev_mx != 0 && prev_my != 0 {
            curr_vel.0 = x;
            curr_vel.1 = y;
            let x_dis = x as f32 - prev_mx as f32;
            let y_dis = y as f32 - prev_my as f32;
            curr_vel.2 = if x_dis == 0.0 { curr_vel.2 } else {x_dis*5.0};
            curr_vel.3 = if y_dis == 0.0 { curr_vel.3 } else {y_dis*5.0};
            vx_grid[curr_vel.0 + curr_vel.1*TEX_SIZE] = curr_vel.2;
            vy_grid[curr_vel.0 + curr_vel.1*TEX_SIZE] = curr_vel.3;
          }
          prev_mx = x;
          prev_my = y;
        }
        _ => ()
      }
    }

    // Process fluids
    vx_grid[curr_vel.0 + curr_vel.1*TEX_SIZE] = curr_vel.2;
    vy_grid[curr_vel.0 + curr_vel.1*TEX_SIZE] = curr_vel.3;
    fluid::step_fluid(&mut tex_data[..], &mut vx_grid[..], &mut vy_grid[..], TEX_SIZE as u32, 0.16, 0.0, true);

    // Re buffer texture
    use std::borrow::Cow;
    let raw_tex_2d = glium::texture::RawImage2d {
      data: Cow::from(tex_data.clone()),
      width: TEX_SIZE as u32, height: TEX_SIZE as u32,
      format: glium::texture::ClientFormat::F32, };
    let raw_tex_2d1 = glium::texture::RawImage2d {
      data: Cow::from(vx_grid.clone()),
      width: TEX_SIZE as u32, height: TEX_SIZE as u32,
      format: glium::texture::ClientFormat::F32, };
    let raw_tex_2d2 = glium::texture::RawImage2d {
      data: Cow::from(vy_grid.clone()),
      width: TEX_SIZE as u32, height: TEX_SIZE as u32,
      format: glium::texture::ClientFormat::F32, };
    let texture = glium::texture::Texture2d::new(&display, raw_tex_2d).unwrap();
    let texture1 = glium::texture::Texture2d::new(&display, raw_tex_2d1).unwrap();
    let texture2 = glium::texture::Texture2d::new(&display, raw_tex_2d2).unwrap();

    // Load texture into uniforms
    let uniforms = uniform! { 
      tex: texture.sampled()
        .wrap_function(glium::uniforms::SamplerWrapFunction::Clamp)
        .magnify_filter(glium::uniforms::MagnifySamplerFilter::Linear),
        tex1: texture1.sampled()
          .wrap_function(glium::uniforms::SamplerWrapFunction::Clamp)
          .magnify_filter(glium::uniforms::MagnifySamplerFilter::Linear),
          tex2: texture2.sampled()
            .wrap_function(glium::uniforms::SamplerWrapFunction::Clamp)
            .magnify_filter(glium::uniforms::MagnifySamplerFilter::Linear)
    };

    // Draw
    use glium::Surface;
    let mut target = display.draw();
    target.clear_color(0.0, 0.0, 1.0, 1.0);
    target.draw(&vbo, &indices, &shader, &uniforms, &Default::default()).unwrap();
    target.finish().unwrap();
  }
}
