/*
* SNP_VIS 
* Benjamin VanderSluis
* bvander@cs.umn.edu
* 110818
*/

//import processing.opengl.*;
import fullscreen.*;


final float dt = 0.005;
final float dt2 = dt*dt;
final float font_height = 16;
final float font_width = 8;
//final float NUM_SNIPS = 966976;
final String GENOME_FILE = "/home/slice/sketchbook/snp_vis/large.txt";
static final char[] charset = {'a', 't', 'g', 'c', 'A', 'T', 'G', 'C', 'X', 'Y'};
final int WIDTH = 1600;
final int HEIGHT = 1200;
  //size(screen.width, screen.height, OPENGL);
  //size(400, 400, OPENGL);



Snip snp_ptr;
Genome my_genome;
SnipSphere my_sphere1, my_sphere2, my_sphere3;
SnipHelix my_helix01, my_helix02;
FullScreen fs;

PFont font;

class Snip
{
  //float x, y, vx, vy, ax, ay;
  float theta, phi, r, v_theta, v_phi, v_r;
  boolean off_screen;
  String N;     // nucleotide [ATCG][ATCG] pair
  String ID;  // snp id
  String CHR; // chromosome
  
  Snip(String id, String chr, String n)
  {
     off_screen = false;
     N = n;
     ID = id;
     CHR = chr;
  }
}

class Genome
{
  // VARIABLES
  Stack snps;
  int  snp_count;
  
  Genome()  // CONSTRUCTOR
  {
    snps = new Stack();
    snp_count = 0;
  }
  
  void Load(String genome_file)
  {
    BufferedReader reader = createReader(genome_file);
 
    try {
      String line;
      while ((line = reader.readLine()) != null) {
        String[] words = split(line, '\t');
        snps.add(new Snip(words[0], words[1], words[3]));
        snp_count += 1;
      }
    }
    catch (Exception e) {
      e.printStackTrace(); 
    }
    
     println("loaded " + my_genome.snp_count + " SNPs from " + genome_file);
  }
  
}

class SnipSphere
{
  // VARIABLES
  
  Stack snp_list;
  int pos_x, pos_y, pos_z, pos_radius;
  

  SnipSphere(Genome the_list, int snp_count, int x, int y, int z, int radius){ // CONSTRUCTOR
    snp_list = new Stack();
    //snp_list = the_list;
    pos_x = x; 
    pos_y = y; 
    pos_z = z;
    pos_radius = radius;

	 while(!the_list.snps.empty() && snp_count > 0){
		snp_ptr = (Snip)the_list.snps.pop();

      snp_ptr.theta = random(0, PI); 
      snp_ptr.phi = random(0,TWO_PI); 
      snp_ptr.r = pos_radius;

      snp_ptr.v_theta = random(0, PI);
      snp_ptr.v_phi = random(0, PI);
      snp_ptr.v_r = 0; // Stay on the surface for now

		snp_list.push(snp_ptr);
		snp_count--;
    }
  }
  
  void update(){
    ListIterator vitr = snp_list.listIterator();
    while(vitr.hasNext()){
      snp_ptr = (Snip)vitr.next();

      snp_ptr.theta = (snp_ptr.theta + snp_ptr.v_theta*dt)%(TWO_PI);
      snp_ptr.phi   = (snp_ptr.phi   + snp_ptr.v_phi  *dt)%(TWO_PI);
      snp_ptr.r     += snp_ptr.v_r    * dt;
    }
  }

  void draw(){
    pushMatrix();
    translate(pos_x, pos_y, pos_z);
    ListIterator vitr = snp_list.listIterator();
    while(vitr.hasNext()){
      snp_ptr = (Snip)vitr.next();
      float x = cos(snp_ptr.theta) * sin(snp_ptr.phi) * snp_ptr.r;
      float y = sin(snp_ptr.theta) * sin(snp_ptr.phi) * snp_ptr.r;
      float z = cos(snp_ptr.phi)   * snp_ptr.r;
  
       //pushMatrix();
       //translate(x, y, z);
       //text(snp_ptr.N, 0, 0, 0);
       //popMatrix();
       text(snp_ptr.N, x, y, z);

    }
    popMatrix();
  }
}

class SnipHelix{
  // VARIABLES
  Stack snp_list;
  int prime_x, prime_y, prime_z;
  int prime_xx, prime_yy, prime_zz;
  int radius;
  float d_theta, the_length;

  SnipHelix(Genome the_list, int snp_count, int x, int y, int z, int xx, int yy, int zz, int rr, float ss){ // CONSTRUCTOR
    snp_list = new Stack();
    prime_x = x;
    prime_y = y;
    prime_z = z;
    prime_xx = xx;
    prime_yy = yy;
    prime_zz = zz;
    radius = rr;
    d_theta = ss;

    int i=0;
    int final_count = snp_count;
    while(!the_list.snps.empty() && snp_count > 0){
      snp_ptr = (Snip)the_list.snps.pop();
      snp_ptr.theta = (float)i/final_count; // [0 - 1) position on strand
      //println("theta "+snp_ptr.theta+"\n");
		snp_list.push(snp_ptr);
		snp_count--;
		i++;
    }

    the_length = sqrt(pow(prime_xx-prime_x,2)+pow(prime_yy-prime_y,2)+pow(prime_zz-prime_z,2));
  }


  void update(){
    ListIterator vitr = snp_list.listIterator();
    while(vitr.hasNext()){
        snp_ptr = (Snip)vitr.next();
        snp_ptr.theta = (snp_ptr.theta + d_theta*dt) % 1;
    }
  }

  void draw(){
    pushMatrix();
    translate(prime_xx, prime_yy, prime_zz);
    if(prime_x - prime_xx != 0){
	    rotateZ(-atan((prime_y - prime_yy) / (prime_x - prime_xx)));
	 }
    if(prime_z - prime_zz != 0){
       rotateX(-atan((prime_y - prime_yy) / (prime_z - prime_zz)));
    }
		// now the helix is vertical from the origin?

    ListIterator vitr = snp_list.listIterator();
    while(vitr.hasNext()){
      snp_ptr = (Snip)vitr.next();
		float num_turns = 5;
      float x = cos(snp_ptr.theta * TWO_PI * num_turns)*radius;
      float y = the_length*(snp_ptr.theta);
      float z = sin(snp_ptr.theta * TWO_PI * num_turns)*radius;

      pushMatrix();
      translate(x, y, z);
      text(snp_ptr.N, 0, 0, 0);
      popMatrix();
    }

    popMatrix();
  }

}

void setup(){


  size(WIDTH, HEIGHT, P3D);
  fs = new FullScreen(this);
  fs.enter();

  my_genome = new Genome();
  my_genome.Load(GENOME_FILE);

  my_sphere1 = new SnipSphere(my_genome, 20, WIDTH*3/4, HEIGHT*3/4, -1000, 1000);
  my_sphere2 = new SnipSphere(my_genome, 40, WIDTH*1/4, HEIGHT*1/4, 0, 100);
  my_sphere3 = new SnipSphere(my_genome, 40, WIDTH*1/8, HEIGHT*4/5, -1000, 100);

  my_helix01 = new SnipHelix(my_genome, 33, (int)WIDTH/8, (int)WIDTH/8, 0, (int)WIDTH/4, (int)(WIDTH/2.28), -(int)WIDTH/2,  30, 0.2);
  my_helix02 = new SnipHelix(my_genome, 33, (int)WIDTH*7/8, 0, (int)-WIDTH/8, (int)WIDTH*5/8, (int)WIDTH/4, 0,  30, 0.4);
  
  
  hint(ENABLE_NATIVE_FONTS);
  textFont(createFont("Sans", 16, false, charset));
  
  background(0);

}


void draw(){
  background(0);
  fill(050,222,050);
  my_sphere1.draw();
  my_sphere1.update();
  my_sphere2.draw();
  my_sphere2.update();
  my_sphere3.draw();
  my_sphere3.update();
  my_helix01.draw();
  my_helix01.update();
  my_helix02.draw();
  my_helix02.update();
}


