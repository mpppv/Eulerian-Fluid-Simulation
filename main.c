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

#define WINDOW_WIDTH        600
#define WINDOW_HEIGHT       900

#define SPHERE_RAD          50

#define CELL_SIZE           5
#define DRAW_GRID           false
#define LINE_WIDTH          1
#define GRID_COLOR          0x262626

#define CELL_SIZE_M         1 //((double)CELL_SIZE / (double)WINDOW_WIDTH)
#define TIME_STEP_HZ        60.0
#define GRAVITY             9.81
#define DENSITY             (double)1000

#define CONVERGENCE_N       1
#define OVERRELAXATION      1

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
double interpolate(double x, double y, char field);
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
double cell_size_m = CELL_SIZE_M;//(double)CELL_SIZE / (double)WINDOW_WIDTH;
double cell_per_s = TIME_STEP_HZ * (double)CELL_SIZE / (double)WINDOW_WIDTH;
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
            else if(event.button.button == SDL_BUTTON_LEFT){//SDL_GetModState() & KMOD_SHIFT
                SDL_GetMouseState(&mouse_x, &mouse_y);
                update_obj(mouse_x/CELL_SIZE, rows - mouse_y/CELL_SIZE);
            }
        }
        
        add_obj_to_cells();
        
        add_gravity();
        for(int i = 0; i < CONVERGENCE_N; ++i){
            projection();
        }
        extrapolate();
        for(int i = 0; i < CONVERGENCE_N; ++i){
            advection();
        }
//        for(int i = rows-1; i >= 0; --i){
//        for(int j = 0; j < columns; ++j){
//        printf("%.2f\t", cells[i][j].v);
//        }printf("\t\t");
//        for(int j = 0; j < columns; ++j){
//        printf("%.2f\t", cells[i][j].u);
//        }printf("\n"); }printf("\n");
//        printf("-----------------------------------------------------------------------------------");
//        printf("-----------------------------------------------------------------------------------\n");
        SDL_FillRect(surface, &black_rect, COLOR_BLACK);
        draw_cells(surface);
        //if(DRAW_GRID) draw_grid(surface);
        draw_sphere(surface);
        
        SDL_UpdateWindowSurface(window);
        SDL_Delay(16.67);
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
            obj.ocup_cll[i][j] = (i - rows/2.0)*(i - rows/2.0) + (j - columns/2.0)*(j - columns/2.0) <= radius * radius;
    }
}

void set_cell_to(cell* cell, bool type, double v, double u){
    cell->is_liquid = type;
    cell->v = v;
    cell->u = u;
    cell->p = 0;
}

void update_obj(int x, int y){
    y = (y < radius + 2)? radius+2 : ((y >= rows    - radius-2)? rows    - radius - 3 : y);
    x = (x < radius + 2)? radius+2 : ((x >= columns - radius-2)? columns - radius - 3 : x);
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
                set_cell_to(&(cells[i][j]), true, (double)(obj.y - cells[i][j].y) * cell_per_s, (double)(obj.x - cells[i][j].x) * cell_per_s);
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
    int s1, s2, s3, s4;
    double s;
    for(int i = 1; i < rows-1; ++i) {
        for(int j = 1; j < columns-1; ++j) {
            if(!cells[i][j].is_liquid) continue;
            
            d = cells[i][j+1].u - cells[i][j].u + cells[i+1][j].v - cells[i][j].v;
            d *= OVERRELAXATION;
            
            s1 = cells[i-1][j].is_liquid? 1 : 0;
            s2 = cells[i+1][j].is_liquid? 1 : 0;
            s3 = cells[i][j-1].is_liquid? 1 : 0;
            s4 = cells[i][j+1].is_liquid? 1 : 0;
            s = s1 + s2 + s3 + s4;
            if(s == 0) continue;
            
            if(s1 == 1) cells[i][j].v   -= d / s;
            if(s2 == 1) cells[i+1][j].v += d / s;
            if(s3 == 1) cells[i][j].u   -= d / s;
            if(s4 == 1) cells[i][j+1].u += d / s;
            //printf("%---.2f\t", d);
            //printf("(%.2f, %.2f) ", cells[i][j].u, cells[i][j].v);
            
            cells[i][j].p += d * ds / s;
        }
        //printf("\n");
    }
}

void extrapolate(void) {
    for (int i = 0; i < columns; ++i) {
        cells[0][i].v = cells[1][i].v;
        cells[rows - 1][i].v = cells[rows - 2][i].v;
    }
    for(int j = 0; j < rows; ++j){
        cells[j][0].u = cells[j][1].u;
        cells[j][columns - 1].u = cells[j][columns - 2].u;
    }
}

double interpolate(double x, double y, char field) {
    
    int x0 = floor(x);
    int y0 = floor(y);
    
    double tx = (x - (double)x0) / cell_size_m;
    double ty = (y - (double)y0) / cell_size_m;
    //printf("%.2f,\t%.2f\n", (x - (double)x0) / ((double)CELL_SIZE / (double)WINDOW_WIDTH), ty);
    //printf("%.2f,\t%.2f\n", x, y);
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
    double x, y;
    
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j)
            copy[i][j] = cells[i][j];
    
    
    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < columns; ++j){
            if(cells[i][j].is_liquid && j < columns - 1){
                v_bar = (copy[i][j].v + copy[i][j+1].v + copy[i-1][j].v + copy[i-1][j+1].v) * 0.25;
                
                x = j * cell_size_m - copy[i][j].u / TIME_STEP_HZ;
                y = (i) * cell_size_m - v_bar / TIME_STEP_HZ;
                y = (y < 0)? 0 : ((y >= rows)?      rows    - 1 : y);
                x = (x < 0)? 0 : ((x >= columns)?   columns - 1 : x);
                
                cells[i][j].u = interpolate(x, y, 'U');
            }
            
            if(cells[i][j].is_liquid && i < rows - 1){
                u_bar = (copy[i][j].u + copy[i+1][j].u + copy[i][j-1].u + copy[i+1][j-1].u) * 0.25;
                
                x = (j) * cell_size_m - u_bar/TIME_STEP_HZ;
                y = i * cell_size_m -  copy[i][j].v/TIME_STEP_HZ;
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
    double conts = 2 / 1.1339745962;
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

    for(int i = rows - 2; i >= 1; --i){
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
                velocity = 0xff;//* (vel - min) / d;
                color = velocity * conts;//(velocity << 16) + (velocity << 8) + velocity;
                
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
