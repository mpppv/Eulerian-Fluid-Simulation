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
#define COLOR_GREY          0x0f0f0f0f
#define LIQUID_COLOR        0xa2c0f5
#define SOLID_COLOR         0xbea1bf

#define PRESSURE_MODE       0

#define WINDOW_WIDTH        900
#define WINDOW_HEIGHT       900

#define CELL_SIZE           5
#define DRAW_GRID           false
#define LINE_WIDTH          2
#define GRID_COLOR          0x262626

#define CELL_SIZE_M         1
#define TIME_STEP_HZ        100
#define GRAVITY             9.81
#define DENSITY             1000

#define CONVERGENCE_N       100
#define OVERRELAXATION      1.9

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
void projection(void);
void advection(void);

void draw_grid(SDL_Surface* surface);
void draw_cell(SDL_Surface* surface, int x, int y, Uint32 color);
void draw_cells(SDL_Surface* surface);


SDL_Rect black_rect = (SDL_Rect){0, 0, WINDOW_WIDTH, WINDOW_HEIGHT};
solid_obj obj;
cell cells[WINDOW_HEIGHT / CELL_SIZE][WINDOW_WIDTH / CELL_SIZE];
cell copy[WINDOW_HEIGHT / CELL_SIZE][WINDOW_WIDTH / CELL_SIZE];

int rows = WINDOW_HEIGHT / CELL_SIZE;
int columns = WINDOW_WIDTH / CELL_SIZE;
int radius = 60 / CELL_SIZE;
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
            else if(event.button.button == SDL_BUTTON_LEFT){ //SDL_GetModState() & KMOD_SHIFT
                SDL_GetMouseState(&mouse_x, &mouse_y);
                if(SDL_GetModState() & KMOD_SHIFT){
                    mouse_y = rows - mouse_y/CELL_SIZE;
                    mouse_x = mouse_x/CELL_SIZE;
                    mouse_y = (mouse_y < radius)? radius : ((mouse_y >= rows - radius)? rows - radius - 1 : mouse_y);
                    for (int i = mouse_y - radius; i < mouse_y + radius + 1; ++i) 
                        cells[i][mouse_x].v += 10000000;
                }
                else
                    update_obj(mouse_x/CELL_SIZE, rows - mouse_y/CELL_SIZE);
            }
        }
        
        SDL_FillRect(surface, &black_rect, COLOR_BLACK);
        add_obj_to_cells();

        //add_gravity();
        for(int i = 0; i < CONVERGENCE_N; ++i)
            projection();
        //for(int i = 0; i < CONVERGENCE_N; ++i)
            advection();
        
        
        draw_cells(surface);
        //if(DRAW_GRID) draw_grid(surface);
        
        SDL_UpdateWindowSurface(window);
        SDL_Delay(10);
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
    y = (y < radius)? radius : ((y >= rows    - radius)? rows    - radius - 1 : y);
    x = (x < radius)? radius : ((x >= columns - radius)? columns - radius - 1 : x);
    obj.v = (y - obj.y) * cell_per_s;
    obj.u = (x - obj.x) * cell_per_s;
    obj.y = y;
    obj.x = x;
    //if(!obj.ocup_cll[y][x]) return;

    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < columns; ++j)
            obj.ocup_cll[i][j] = false;

    for (int i = y - radius; i < y + radius + 1; ++i) 
        for (int j = x - radius; j < x + radius + 1; ++j)
            obj.ocup_cll[i][j] = (i - y)*(i - y) + (j - x)*(j - x) <= radius * radius;
}

void add_obj_to_cells(void){
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < columns; ++j) {
            if(obj.ocup_cll[i][j] && cells[i][j].is_liquid) 
                set_cell_to(&(cells[i][j]), false, obj.v, obj.u);
            else if(!obj.ocup_cll[i][j] && !cells[i][j].is_liquid)
                set_cell_to(&(cells[i][j]), true, (obj.y - cells[i][j].y) * cell_per_s, (obj.x - cells[i][j].x) * cell_per_s); ///TODO
        }
}



void add_gravity(void){
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < columns; ++j)
            if(cells[i][j].is_liquid)
                cells[i][j].v += g_accel_step;
}

void projection(void){
    double d;
    double ds = DENSITY * cell_per_s;
    int s, s1, s2, s3, s4;
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < columns; ++j)
            if(cells[i][j].is_liquid){
                if(i == rows-1 && j == columns-1)
                    d = -cells[i][j].u - cells[i][j].v;
                else if(i == rows-1)
                    d = cells[i][j+1].u - cells[i][j].u - cells[i][j].v;
                else if(j == columns-1)
                    d = -cells[i][j].u + cells[i+1][j].v - cells[i][j].v;
                else
                    d = cells[i][j+1].u - cells[i][j].u + cells[i+1][j].v - cells[i][j].v;
                d *= OVERRELAXATION;
                
                s1 = (i > 0)? 1 : 0            ;//&& cells[i-1][j].is_liquid;
                s2 = (i < rows-1) ? 1 : 0       ;//&& cells[i+1][j].is_liquid;
                s3 = (j > 0)  ? 1 : 0           ;//&& cells[i][j-1].is_liquid;
                s4 = (j < columns-1) ? 1 : 0    ;//&& cells[i][j+1].is_liquid;
                s = s1 + s2 + s3 + s4;
                if(s) continue;
                
                if(s1 == 1) cells[i][j].v -= d / s;
                if((s2 == 1) && cells[i+1][j].is_liquid) cells[i+1][j].v += d / s;
                if(s3 == 1) cells[i][j].u -= d / s;
                if((s4 == 1) && cells[i][j+1].is_liquid) cells[i][j+1].u += d / s;
                
                cells[i][j].p += d * ds / s;
            }
}

void advection(void){
    double v_bar;
    int x, y;
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j)
            copy[i][j] = cells[i][j];
    }
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < columns; ++j){
            if(!copy[i][j].is_liquid) continue;
//            
//            if(i == 0 && j == columns-1)
//                v_bar = cells[i][j].v;
//            else if(i == 0)
//                v_bar = (cells[i][j].v + cells[i][j+1].v) / 2;
//            else if(j == columns-1)
//                v_bar = (cells[i][j].v + cells[i-1][j].v) / 2;
//            else
//                v_bar = (cells[i][j].v + cells[i][j+1].v + cells[i-1][j].v + cells[i-1][j+1].v) / 4;
//            y = i + v_bar / ms;
//            x = j + copy[i][j].u / ms;
//            y = (y < 0)? 0 : ((y >= rows)?      rows    - 1 : y);
//            x = (x < 0)? 0 : ((x >= columns)?   columns - 1 : x);
//            
//            if(cells[i][j].is_liquid) cells[i][j] = copy[y][x];
            
             y = i + copy[i][j].v / cell_per_s;
             x = j + copy[i][j].u / cell_per_s;
             y = (y < 0)? 0 : ((y >= rows)? rows - 1 : y);
             x = (x < 0)? 0 : ((x >= columns)? columns - 1 : x);
            
             if(cells[y][x].is_liquid) cells[y][x] = copy[i][j];
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
                if(cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u > max) max = cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u;
                else if(cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u < min) min = cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u;
            }
        }
        max = sqrt(max); min = sqrt(min);
    }
    double d = max - min;

    for(int i = 1; i < rows-1; ++i){
        for(int j = 1; j < columns-1; ++j){
            if(cells[i][j].is_liquid){
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
                    velocity = 0xff * (vel - min) / d;
                    color = velocity * 2 / 1.1339745962;
                    //green
                    if(cells[i][j].u >= 0 && cos_v >= 0){
                        color *= cos_v;
                        if(color <= velocity)    color = (velocity << 16) + (color << 8);
                        else                color = (velocity << 16)+(velocity << 8) - ((color - velocity) << 16);
                    }
                    else if(cells[i][j].u < 0 && cos_v >= 0.8660254038){ // √3/2
                        color *= (2 - cos_v);
                        color = (velocity << 16)+(velocity << 8) - ((color - velocity) << 16);
                    }
                    //blue
                    else if(cells[i][j].u < 0 && cos_v >= -0.8660254038){
                        color = (0.8660254038 - cos_v) * velocity * 2 / 1.7320508076;
                        if(color <= velocity)    color = (velocity << 8) + (color);
                        else                color = (velocity << 8)+(velocity) - ((color - velocity) << 8);
                    }
                    //red
                    else if(cells[i][j].u < 0){
                        color *= - (0.8660254038 + cos_v);
                        color = velocity + (color << 16);
                    }
                    else{
                        color *= (1.1339745962 + cos_v);
                        if(color <= velocity)    color = velocity + (color << 16);
                        else                color = (velocity << 16)+(velocity) - (color - velocity);
                    }
                    draw_cell(surface, j, i, color);
                }
            }
            else draw_cell(surface, j, i, SOLID_COLOR);
        }
    }
}
