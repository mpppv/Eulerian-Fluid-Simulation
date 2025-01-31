//
//  main.c
//  Fluid_Simulation
//
//  Created by mpppv on 2025/01/26.
//
#include <SDL2/SDL.h>
#include <math.h>
#include <stdbool.h>

#define COLOR_BLACK         0x00000000
#define SOLID_COLOR         0xbea1bf

#define PRESSURE_MODE       false

#define WINDOW_WIDTH        900
#define WINDOW_HEIGHT       900

#define SPHERE_RAD          50

#define CELL_SIZE           5
#define DRAW_GRID           false
#define LINE_WIDTH          2
#define GRID_COLOR          0x262626

#define CELL_SIZE_M         1 //(CELL_SIZE / WINDOW_WIDTH)
#define TIME_STEP_HZ        10
#define GRAVITY             9.81
#define DENSITY             1000

#define CONVERGENCE_N       100
#define OVERRELAXATION      1.8

typedef struct{
    bool is_liquid;
    int x, y;
    
    double v, u;
    double p;
}cell;

typedef struct{
    int x, y;
    double v, u;
    bool ocup_cll[WINDOW_HEIGHT / CELL_SIZE][WINDOW_WIDTH / CELL_SIZE];
}solid_obj;

void init_cells(void);
void init_obj(void);
void set_cell_to(cell* cell, bool type, double v, double u);
void update_obj(int x, int y);

void add_obj_to_cells(void);

void add_gravity(void);
void extrapolate(void);
double interpolate(int x, int y, char field);
void projection(void);
void advection(void);

void draw_grid(SDL_Surface* surface);
void draw_cell(SDL_Surface* surface, int x, int y, Uint32 color);
void draw_cells(SDL_Surface* surface);
void draw_sphere(SDL_Surface* surface);


SDL_Rect black_rect = (SDL_Rect){0, 0, WINDOW_WIDTH, WINDOW_HEIGHT};
solid_obj obj;
cell cells[WINDOW_HEIGHT / CELL_SIZE][WINDOW_WIDTH / CELL_SIZE];
cell copy[WINDOW_HEIGHT / CELL_SIZE][WINDOW_WIDTH / CELL_SIZE];

int rows = WINDOW_HEIGHT / CELL_SIZE;
int columns = WINDOW_WIDTH / CELL_SIZE;
int radius = SPHERE_RAD / CELL_SIZE;
double cell_per_s = CELL_SIZE_M * TIME_STEP_HZ;
double g_accel_step = GRAVITY / TIME_STEP_HZ;


int main(void) {
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Init(SDL_INIT_EVENTS);
    SDL_Window* window = SDL_CreateWindow("Liquid Simulation",
                                          SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                          WINDOW_WIDTH,           WINDOW_HEIGHT, SDL_WINDOW_METAL);

    SDL_Surface* surface = SDL_GetWindowSurface(window);

    init_obj();
    init_cells();
    
    bool running = true;
    int mouse_x = 0, mouse_y = 0;
    SDL_Event event;
    while(running){
        while(SDL_PollEvent(&event)){
            if(event.type == SDL_QUIT) running = false;
            else if(event.button.button == SDL_BUTTON_LEFT){
                SDL_GetMouseState(&mouse_x, &mouse_y);
                if(SDL_GetModState() & KMOD_SHIFT){
                    mouse_y = rows - mouse_y/CELL_SIZE;
                    mouse_x = mouse_x/CELL_SIZE;
                    mouse_y = (mouse_y < radius)? radius : ((mouse_y >= rows - radius)? rows - radius - 1 : mouse_y);
                    for (int i = mouse_y - radius; i < mouse_y + radius + 1; ++i) 
                        cells[i][mouse_x].u += 10;
                }
                else
                    update_obj(mouse_x/CELL_SIZE, rows - mouse_y/CELL_SIZE);
            }
        }
        
        add_obj_to_cells();

        add_gravity();
        for(int i = 0; i < CONVERGENCE_N; ++i)
            projection();
        extrapolate();
        advection();
        
        SDL_FillRect(surface, &black_rect, COLOR_BLACK);
        draw_cells(surface);
        //if(DRAW_GRID) draw_grid(surface);
        draw_sphere(surface);
        
        SDL_UpdateWindowSurface(window);
        SDL_Delay(17);
    }
}


void init_cells(void){
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j)
            cells[i][j] = (cell){true, j, i, 0, 0, 0};
}

void init_obj(void){
    obj.x = columns/2; obj.y = rows/2; obj.v = obj.u = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j)
            obj.ocup_cll[i][j] = (i - rows/2)*(i - rows/2) + (j - columns/2)*(j - columns/2) <= radius * radius;
    }
}

void set_cell_to(cell* cell, bool type, double v, double u){
    cell->is_liquid = type;
    cell->v = v;
    cell->u = u;
    cell->p = 0;
}

void update_obj(int x, int y){
    y = (y < radius + 1)? radius+1 : ((y >= rows    - radius-1)? rows    - radius - 2 : y);
    x = (x < radius + 1)? radius+1 : ((x >= columns - radius-1)? columns - radius - 2 : x);
    obj.v = (y - obj.y) * cell_per_s;
    obj.u = (x - obj.x) * cell_per_s;
    obj.y = y;
    obj.x = x;

    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < columns; ++j)
            obj.ocup_cll[i][j] = false;

    for (int i = y - radius; i < y + radius + 1; ++i) 
        for (int j = x - radius; j < x + radius + 1; ++j)
            obj.ocup_cll[i][j] = (i - y)*(i - y) + (j - x)*(j - x) <= radius * radius;
}

void add_obj_to_cells(void){
    for(int i = 1; i < rows-1; ++i)
        for(int j = 1; j < columns-1; ++j) {
            if(obj.ocup_cll[i][j] && cells[i][j].is_liquid)
                set_cell_to(&(cells[i][j]), false, obj.v, obj.u);
            else if(!obj.ocup_cll[i][j] && !cells[i][j].is_liquid)
                set_cell_to(&(cells[i][j]), true, obj.v, obj.u);//(obj.y - cells[i][j].y) * cell_per_s, (obj.x - cells[i][j].x) * cell_per_s);
        }
}



void add_gravity(void){
    for(int i = 1; i < rows-1; ++i)
        for(int j = 1; j < columns-1; ++j)
            if(cells[i][j].is_liquid && cells[i - 1][j].is_liquid)
                cells[i][j].v -= g_accel_step;
}

void projection(void){
    double d;
    double ds = DENSITY * cell_per_s;
    int s, s1, s2, s3, s4;
    for(int i = 1; i < rows-1; ++i) {
        for(int j = 1; j < columns-1; ++j) {
            if(cells[i][j].is_liquid){
                d = cells[i][j+1].u - cells[i][j].u + cells[i+1][j].v - cells[i][j].v;
                d *= OVERRELAXATION;
                
                s1 = (i > 0              && cells[i-1][j].is_liquid)? 1 : 0;
                s2 = (i < rows-1         && cells[i+1][j].is_liquid)? 1 : 0;
                s3 = (j > 0              && cells[i][j-1].is_liquid)? 1 : 0;
                s4 = (j < columns-1      && cells[i][j+1].is_liquid)? 1 : 0;
                s = s1 + s2 + s3 + s4;
                if(s) continue;
                
                if(s1 == 1) cells[i][j].v -= d / s;
                if(s2 == 1/* && cells[i+1][j].is_liquid*/) cells[i+1][j].v += d / s;
                if(s3 == 1) cells[i][j].u -= d / s;
                if(s4 == 1/* && cells[i][j+1].is_liquid*/) cells[i][j+1].u += d / s;
                
                cells[i][j].p += d * ds / s;
            }
        }
    }
}

void extrapolate(void) {
    for (int i = 0; i < columns; ++i) {
        cells[0][i].u = cells[1][i].u;
        cells[rows - 1][i].u = cells[rows - 2][i].u;
    }
    for(int j = 0; j < rows; ++j){
        cells[j][0].v = cells[j][1].v;
        cells[j][columns - 1] = cells[j][columns - 2];
    }
}

double interpolate(int x, int y, char field) {
    double dx = 0, dy = 0;
    switch (field) {
        case 'V': dx = CELL_SIZE_M / 2; break;
        case 'U': dy = CELL_SIZE_M / 2; break;
        case 'S': dx = dy =  CELL_SIZE_M / 2; break;
        default: break;
    }
    
    //
    
    int x0 = (x-dx)/CELL_SIZE_M < columns-1 ? (x-dx)/CELL_SIZE_M : columns-1;
    int y0 = (y-dy)/CELL_SIZE_M < rows-1    ? (y-dy)/CELL_SIZE_M : rows-1;
    
    int tx = ((x-dx) - x0*CELL_SIZE_M) / CELL_SIZE_M;
    int ty = ((y-dy) - y0*CELL_SIZE_M) / CELL_SIZE_M;
    
    int x1 = x0 + 1 < columns-1 ? x0 + 1 : columns-1;
    int y1 = y0 + 1 < rows-1    ? y0 + 1 : rows-1;
    
    if(field == 'V')
        return (1-tx)   * (1-ty)    * copy[y0][x0].v +
               tx       * (1-ty)    * copy[y0][x1].v +
               tx       * ty        * copy[y1][x1].v +
               (1-tx)   * ty        * copy[y1][x0].v;
    else if(field == 'U')
        return (1-tx)   * (1-ty)    * copy[y0][x0].u +
               tx       * (1-ty)    * copy[y0][x1].u +
               tx       * ty        * copy[y1][x1].u +
               (1-tx)   * ty        * copy[y1][x0].u;
    else //if(field == 'S')
        return (1-tx)   * (1-ty)    * copy[y0][x0].is_liquid +
               tx       * (1-ty)    * copy[y0][x1].is_liquid +
               tx       * ty        * copy[y1][x1].is_liquid +
               (1-tx)   * ty        * copy[y1][x0].is_liquid;
}

void advection(void){
    double v_bar, u_bar;
    int x, y;
    
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j)
            copy[i][j] = cells[i][j];
    
    
    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < columns; ++j){
            if(cells[i][j].is_liquid && i < rows - 1){
                v_bar = (copy[i][j].v + copy[i][j+1].v + copy[i-1][j].v + copy[i-1][j+1].v) / 4;
                
                x = j * CELL_SIZE_M - copy[i][j].u/TIME_STEP_HZ;
                y = (i + 0.5) * CELL_SIZE_M -  v_bar/TIME_STEP_HZ;
                y = (y < 0)? 0 : ((y >= rows)?      rows    - 1 : y);
                x = (x < 0)? 0 : ((x >= columns)?   columns - 1 : x);
                
                cells[i][j].u = interpolate(x, y, 'U');
            }
            
            if(cells[i][j].is_liquid && j < columns - 1){
                u_bar = (copy[i][j].u + copy[i+1][j].u + copy[i][j-1].u + copy[i+1][j-1].u) / 4;
                
                x = (j + 0.5) * CELL_SIZE_M - u_bar/TIME_STEP_HZ;
                y = i * CELL_SIZE_M -  copy[i][j].v/TIME_STEP_HZ;
                y = (y < 0)? 0 : ((y >= rows)?      rows    - 1 : y);
                x = (x < 0)? 0 : ((x >= columns)?   columns - 1 : x);
                
                cells[i][j].v = interpolate(x, y, 'V');
            }
        }
    }
}



void draw_grid(SDL_Surface* surface){
    for(int i = 0; i < columns; ++i){
        SDL_Rect column = (SDL_Rect){i * CELL_SIZE, 0, LINE_WIDTH, WINDOW_HEIGHT};
        SDL_FillRect(surface, &column, GRID_COLOR);
    }
    for(int i = 0; i < rows; ++i){
        SDL_Rect row = (SDL_Rect){0, i * CELL_SIZE, WINDOW_WIDTH, LINE_WIDTH};
        SDL_FillRect(surface, &row, GRID_COLOR);
    }
}

void draw_cell(SDL_Surface* surface, int x, int y, Uint32 color){
    SDL_Rect cell = (SDL_Rect){(x)*CELL_SIZE, (rows-1-y)*CELL_SIZE, CELL_SIZE, CELL_SIZE};
    SDL_FillRect(surface, &cell, color);
}
void draw_cells(SDL_Surface* surface){
    Uint32 color;
    int pressure, velocity;
    //int x_max = 0, y_max = 0, x_min = 0, y_min = 0;
    double vel, cos_v;
    double min = 1.797693e+308, max = -1.797693e+308;
    if(PRESSURE_MODE){
        for(int i = 1; i < rows-1; ++i){
            for(int j = 1; j < columns-1; ++j){
                if(cells[i][j].p > max) max = cells[i][j].p;
                else if(cells[i][j].p < min) min = cells[i][j].p;
            }
        }
    }
    else{
        for(int i = 1; i < rows-1; ++i){
            for(int j = 1; j < columns-1; ++j){
                if(!cells[i][j].is_liquid) continue;
                if(cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u > max){ max = cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u;} //x_max= j; y_max = i;}
                else if(cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u < min){ min = cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u;} //x_min= j; y_min = i;}
            }
        }
        max = sqrt(max); min = sqrt(min);
    }
    double d = max - min;
    for(int i = 1; i < rows-1; ++i){
        for(int j = 1; j < columns-1; ++j){
            if(PRESSURE_MODE) {
                pressure = 3 * 255 * (cells[i][j].p - min) / d;
                if (pressure >= 510)        color = 0xff00ff - pressure + 510;
                else if (pressure >= 255)   color = 0x0000ff + ((pressure - 255) << 16);
                else                        color = 0x00ffff - (pressure << 8);
                draw_cell(surface, j, i, color);
            }
            else{
                vel = sqrt(cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u);
                cos_v = cells[i][j].v / vel;
                velocity = 0xff * (vel) / 100;
                color = velocity * 2 / 1.1339745962;//(velocity << 16) + (velocity << 8) + velocity;
                
                //green
                if(cells[i][j].u >= 0 && cos_v >= 0){
                    color *= cos_v;
                    if(color <= velocity)   color = (velocity << 16) + (color << 8);
                    else                    color = (velocity << 16)+(velocity << 8) - ((color - velocity) << 16);
                }
                else if(cells[i][j].u < 0 && cos_v >= 0.8660254038){ // √3/2
                    color *= (2 - cos_v);
                    color = (velocity << 16)+(velocity << 8) - ((color - velocity) << 16);
                }
                //blue
                else if(cells[i][j].u < 0 && cos_v >= -0.8660254038){
                    color = (0.8660254038 - cos_v) * velocity * 2 / 1.7320508076;
                    if(color <= velocity)   color = (velocity << 8) + color;
                    else                    color = (velocity << 8)+velocity - ((color - velocity) << 8);
                }
                //red
                else if(cells[i][j].u < 0){
                    color *= - (0.8660254038 + cos_v);
                    color = velocity + (color << 16);
                }
                else{
                    color *= (1.1339745962 + cos_v);
                    if(color <= velocity)   color = velocity + (color << 16);
                    else                    color = (velocity << 16)+velocity - (color - velocity);
                }
                
                draw_cell(surface, j, i, color);
            }
        }
    }
    //printf("%f,\t%f\n", max, min);
    //draw_cell(surface, x_max, y_max, 0xeb5b5b);
    //draw_cell(surface, x_min, y_min, 0x5bd3eb);
}

void draw_sphere(SDL_Surface* surface) {
    double rad_sq = (SPHERE_RAD + CELL_SIZE) * (SPHERE_RAD + CELL_SIZE);
    double rad_sq_obr = (SPHERE_RAD + CELL_SIZE*1.2) * (SPHERE_RAD + CELL_SIZE*1.2);
    double rad_sq_surp = (SPHERE_RAD + CELL_SIZE*2) * (SPHERE_RAD + CELL_SIZE*2);
    double dist;
    int x = obj.x * CELL_SIZE, y = obj.y * CELL_SIZE;
    for(int i = y - SPHERE_RAD - CELL_SIZE*2; i < y + SPHERE_RAD + CELL_SIZE*2; ++i)
        for(int j = x - SPHERE_RAD - CELL_SIZE*2; j < x + SPHERE_RAD + CELL_SIZE*2; ++j){
            dist = (i - y)*(i - y) + (j - x)*(j - x);
            SDL_Rect pix = (SDL_Rect){(j), (WINDOW_HEIGHT-1-i), 1, 1};
            if(dist < rad_sq)
                SDL_FillRect(surface, &pix, 0xf52525);
            else if(dist <= rad_sq_obr)
                SDL_FillRect(surface, &pix, 0x8c1818);
            else if(dist <= rad_sq_surp)
                SDL_FillRect(surface, &pix, 0x2b0707);
        }
}
