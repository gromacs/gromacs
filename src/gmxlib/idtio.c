/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_idtio_c = "$Id$";

/*
 * This module implements synchronous and asynchronous communication via the
 * dual ported ram in the SPC-860 system.
 *
 * Protocol description:
 *
 * Sender and receiver must complete their previous transfer on the specified
 * channel first (busy flag must be cleared for the channel) before they can
 * start a new transfer. A large block of data is split up into a number of
 * blocks which fit into the i/o buffer. A transfer is finished after all
 * blocks have been transferred and the receiver has acknowledged the reception
 * of the last block.
 *
 * Step   Sender                           Receiver
 *
 * 1      gmx_tx(),gmx_txs()               gmx_rx(),gmx_rxs()
 *             .                                .
 * 2      Copy memory to io buffer              .
 *        signal buffer filled                  .
 *             .                                .
 * 3           .                           Wait for buffer filled signal
 *             .                                .
 * 4      Wait for buffer empty signal          .
 *             .                                .
 * 5           .                           Copy io buffer to memory
 *             .                           signal buffer empty
 *             .                                .
 * 6      Buffer transfer done                  .
 *             .                                .
 * 7      Repeat step 2 & 4 until          Repeat step 3 & 5 until
 *        complete transfer done           complete transfer done
 *             .                                .
 * 8      return                           return
 *
 * As is clear from the previous, the sender and the receiver need not to
 * enter the communication routines at the same time. The sender as well as
 * the receiver may invoke its routine at every moment. The protocol
 * guarantuees that no data is lost. The only restriction is that the sender
 * and receiver have the same transfer length. A message with length 0 (zero)
 * is also communicated. Although no data transfer is necessary, the protocol
 * is handled in the same way as a normal message and thus causes the channel
 * to be used for synchronisation. The synchronous and asynchronous routines
 * are interchangable between sender and receiver, so a sender may use the
 * asynchronous version while the receiver uses the synchronous version and
 * vice versa. Naturally both may use the same type of communication
 * primitives.
 *
 * Due to the fact that DMA is not available, all data copying is done by
 * the node. Because interrupts aren't used for efficiency reasons,
 * real asynchronous communication is not possible, messages are copied
 * from and to the io buffers while within the communication routines. So
 * an asynchronous communication is implemented by doing the nonblocking
 * part in the setup and the blocking part in the gmx_wait() routines.
 * Notice that while another asynchronous transfer is busy on the specified
 * channel, starting a new one implies waiting for the previous to finish.
 *
 * IDT layout:
 *
 *             left idt                       right idt
 *
 *   +----------+                 +------+                 +----------+
 *   | 3ff 1023 |  irq (self)     | I860 |    irq (other)  | 7ff 2047 |
 *   +----------+                 |      |                 +----------+
 *   | 3fe 1022 |  irq (other)    |      |    irq (self)   | 7fe 2046 |
 *  -+----------+--------------   |      |   --------------+----------+-
 *   | 3fd 1021 |  rx cmd         |      |    tx cmd       | 7fd 2045 |
 *   +----------+                 |      |                 +----------+
 *   | 3fc 1020 |  tx cmd         |      |    rx cmd       | 7fc 2044 |
 *  -+----------+--------------   |      |   --------------+----------+-------
 *   | 3fb 1019 |                 |      |                 | 7fb 2043 |  ^
 *   +----------+                 |      |                 +----------+  |
 *   .          .                 |      |                 .          .  BSIZE
 *   +----------+                 |      |                 +----------+  |
 *   | 1fe  510 |  rx buffer      |      |    tx buffer    | 5fe 1534 |  v
 *  -+----------+--------------   |      |   --------------+----------+-------
 *   | 1fd  509 |                 |      |                 | 5fd 1533 |  ^
 *   +----------+                 |      |                 +----------+  |
 *   .          .                 |      |                 .          .  BSIZE
 *   +----------+                 |      |                 +----------+  |
 *   | 000    0 |  tx buffer      |      |    rx buffer    | 400 1024 |  v
 *   +----------+                 +------+                 +----------+-------
 *
 * Numbers are offsets relative to IDT base address, names like tx and rx are
 * local for this i860. The actual buffer size may be smaller (aligned for
 * word size, 4 bytes) for speed, see also module inout.s.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <rdklib.h>
#include "fatal.h"
#include "delay.h"
#include "network.h"
#include "inout.h"
#include "idtio.h"
#include "buffer.h"
#include "synclib.h"
#include "main.h"

#define SERVER_VERSION	"SERVER_VERSION"
#define USE_VERSION1	"S860.1"

#define	IDTNR		2	/* Number of available idt's.                */
#define	IDTLEFT		0x0	/* Offset of my idt.                         */
#define	IDTRIGHT	0x400	/* Offset of neighbour nodes idt.       */
#define	IDTSIZE		0x400	/* Total size of idt.                        */
#define	IRQBASE		0x3fe	/* Interrupt locations at the top of an idt  */
#define	BSIZE		510	/* Usable buffer size, see comment above     */
#define	TX		0	/* Index of sender status.                   */
#define	RX		1	/* Index of receiver status.                 */
#define	TXRXNR		2	/* Number of directions (sender & receiver). */
#define	BUF_DONE	0x55	/* Signals buffer transfer is done.          */
#define	BUF_BUSY	((~(BUF_DONE))&0xff)	/* Signals buffer is busy.   */
#define	LINE_WIDTH	16	/* Hex dump bytes per line.                  */
#define MIN(a,b)	(((a)>(b))?(b):(a))
#define WORK_NAME(work)	\
  (((work)==buf_idle)?"buf_idle": \
   ((((work)==buf_accept)?"buf_accept": \
     ((((work)==buf_send)?"buf_send": \
       ((((work)==buf_avail)?"buf_avail":(NULL))))))))

typedef struct t_stat
{
  char *name;	/* Name of the channel (debugging).                          */
  int cmd;	/* Offset of command byte in io ram.                         */
  int buf;	/* Offset of buffer in io ram.                               */
  int len;	/* Actual number of bytes to transfer.                       */
  int bsize;	/* Actual usable buffer size                                 */
  int bcnt;	/* Total block count, part of rx/tx.                         */
  int tcnt;	/* Total transfer count for tx/rx.                           */
  int icnt;	/* Total idle count, incremented while waiting.              */
  char *data;	/* Pointer to data bytes in memory.                          */
  void (*work)(struct t_stat* stat);	/* State of transfer.                */
} t_stat;

static int idt_id=0;
static int idt_busy=0;
static char *stat_names[IDTNR][TXRXNR]={{"txl","rxl"},{"txr","rxr"}};
static t_stat sys_stat[IDTNR][TXRXNR];

static void buf_idle(t_stat *stat);
static void buf_accept(t_stat *stat);
static void buf_send(t_stat *stat);
static void buf_avail(t_stat *stat);

/*
 * Debug routines start
 */
#define	COL_WIDTH	19
#define line(what) \
  do \
    { \
      int i,j,len; \
      for (i=0; i<IDTNR; i++) \
        for (j=0; j<TXRXNR; j++) \
          { \
            len=what; \
            (void) fprintf(fp,"%*s",COL_WIDTH-len,""); \
          } \
      (void) fprintf(fp,"\n"); \
    } \
  while (0)

static void print_connections(FILE *fp,t_stat stat[IDTNR][TXRXNR])
{
  int i,j;

  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      (void) fprintf(fp,"cmd at %s (0x%x) should be connected to %s (0x%x)\n",
                     stat[i][j].name,stat[i][j].cmd,
                     stat[IDTNR-1-i][TXRXNR-1-j].name,
                     stat[IDTNR-1-i][TXRXNR-1-j].cmd);
  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      (void) fprintf(fp,"data at %s (0x%x,%d) should be connected to %s "
                     "(0x%x,%d)\n",
                     stat[i][j].name,stat[i][j].buf,stat[i][j].bsize,
                     stat[IDTNR-1-i][TXRXNR-1-j].name,
                     stat[IDTNR-1-i][TXRXNR-1-j].buf,
                     stat[IDTNR-1-i][TXRXNR-1-j].bsize);
}

static void print_chars(FILE *fp,char s[],int len)
{
  int i;
  
  if (len)
    {
      for (i=len; i<LINE_WIDTH; i++) (void) fprintf(fp,"   ");
      s[len]='\0';
      (void) fprintf(fp," | %s\n",s);
    }
}

static int ascii(int b)
{     
  if ((b<' ')||(b>'~')) return ('.'); else return b;
}

static void idt_dump(FILE *fp,char *title,int offset,int len)
{
  int i,b,index;
  char s[LINE_WIDTH+1];
  
  index=0;
  s[0]='\0';
  if (len) (void) fprintf(fp,"idtdump %s:\n",title);
  for (i=0; i<len; i++)
    {
      b=peek_io_buf(IDT7130,offset);
      if ((index==0)||(index==LINE_WIDTH))
        {
          print_chars(fp,s,index); 
          (void) fprintf(fp,"%.8X :",offset);
          index=0;
        }
      (void) fprintf(fp," %.2X",b);
      s[index++]=ascii(b);
      offset++;
    }
  print_chars(fp,s,index);
}

static int fprintw(FILE *fp,char *title,void (*work)(struct t_stat* stat))
{
  char *work_name;

  if ((work_name=WORK_NAME(work))==NULL)
    return fprintf(fp,"%s0x%x",title,work);
  else
    return fprintf(fp,"%s%s",title,work_name);
}

static void clear_cnt(t_stat stat[IDTNR][TXRXNR])
{
  int i,j;

  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      {
        stat[i][j].bcnt=0;
        stat[i][j].tcnt=0;
        stat[i][j].icnt=0;
      }
}


static void print_status(FILE *fp,int nodeid,char *s,
                         t_stat stat[IDTNR][TXRXNR],t_stat *mark)
{
  (void) fprintf(fp,"idt nodeid: %d, %s\n",nodeid,s);
  line(fprintf(fp,"%s%s:",stat[i][j].name,(mark==&stat[i][j])?" (*)":""));
  line(fprintf(fp,"cmd   = 0x%x",stat[i][j].cmd));
  line(fprintf(fp,"(cmd) = 0x%x",peek_io_buf(IDT7130,stat[i][j].cmd)));
  line(fprintf(fp,"buf   = 0x%x",stat[i][j].buf));
  line(fprintf(fp,"len   = %d",stat[i][j].len));
  line(fprintf(fp,"bsize = %d",stat[i][j].bsize));
  line(fprintf(fp,"bcnt  = %d",stat[i][j].bcnt));
  line(fprintf(fp,"tcnt  = %d",stat[i][j].tcnt));
  line(fprintf(fp,"icnt  = %d",stat[i][j].icnt));
  line(fprintf(fp,"data  = 0x%x",stat[i][j].data));
  line(fprintw(fp,"work  = ",stat[i][j].work));
  fflush(fp);
  clear_cnt(stat);
}

static void init_buf(t_stat *stat,int data)
{
  int i;

  for (i=0; i<stat->bsize; i++) poke_io_buf(IDT7130,stat->buf+i,data);
}

static void init_bufs(t_stat stat[IDTNR][TXRXNR],int data)
{
  int i,j;

  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      init_buf(&stat[i][j],data);
}

static int check_buf(FILE *fp,t_stat *stat,int data)
{
  int i,peeked,err;
  
  err=0;
  for (i=0; i<stat->bsize; i++) 
    if ((peeked=peek_io_buf(IDT7130,stat->buf+i))!=data)
      {
        err++;
        (void) fprintf(fp,"data at 0x%x in %s is 0x%x, should be 0x%x\n",
                       stat->buf+i,stat->name,peeked,data);
      }
  return err;
}

static void init_cmd(t_stat stat[IDTNR][TXRXNR],int cmd)
{
  int i,j;

  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      poke_io_buf(IDT7130,stat[i][j].cmd,cmd);
}  

static void test_connectivity(FILE *fp,t_stat stat[IDTNR][TXRXNR])
{
  int i,j,init,err;

  err=0;
  print_connections(fp,stat);
  (void) fprintf(fp,"short circuit test:\ntesting initialisation ..."); 
  fflush(fp);
  init_cmd(stat,BUF_DONE);
  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      if ((init=peek_io_buf(IDT7130,stat[i][j].cmd))!=BUF_DONE)
        {
          if (err++==0) (void) printf("\n");
          (void) fprintf(fp,"cmd at %s (0x%x) is 0x%x, should be 0x%x\n",
                         stat[i][j].name,stat[i][j].cmd,init,BUF_DONE);
        }
  (void) fprintf(fp," done\ntesting cmd connectivity ..."); 
  fflush(fp);
  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      {
        init_cmd(stat,BUF_DONE); /* All others should have wrong value */
        poke_io_buf(IDT7130,stat[i][j].cmd,BUF_BUSY);
        if ((init=peek_io_buf(IDT7130,stat[IDTNR-1-i][TXRXNR-1-j].cmd))!=
            BUF_BUSY)
          {
            if (err++==0) (void) printf("\n");
            (void) fprintf(fp,"cmd at %s (0x%x) should be connected to"
                           " %s (0x%x)\n",
                           stat[i][j].name,stat[i][j].cmd,
                           stat[IDTNR-1-i][TXRXNR-1-j].name,
                           stat[IDTNR-1-i][TXRXNR-1-j].cmd);
            (void) fprintf(fp,"cmd at %s (0x%x) is 0x%x, should be 0x%x\n",
                           stat[IDTNR-1-i][TXRXNR-1-j].name,
                           stat[IDTNR-1-i][TXRXNR-1-j].cmd,
                           init,BUF_BUSY);
          }
      }
  (void) fprintf(fp," done\ntesting data connectivity ..."); 
  fflush(fp);
  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      {
        init_bufs(stat,0); /* All others should have wrong value */
        init_buf(&stat[i][j],0xff);
        err+=check_buf(fp,&stat[IDTNR-1-i][TXRXNR-1-j],0xff);
        init_bufs(stat,0xff); /* All others should have wrong value */
        init_buf(&stat[i][j],0);
        err+=check_buf(fp,&stat[IDTNR-1-i][TXRXNR-1-j],0);
      }
  (void) fprintf(fp," done\n%d error(s) detected\n",err);
  fflush(fp);
}

#define BUFSIZE	25000

static int com_test(FILE *fp,char *title,int src,int dest)
{
  int buf[2][BUFSIZE];
  
  (void) fprintf(fp,"%s, src=%d, dest=%d\n",title,src,dest);
  clear_buff(buf[0],BUFSIZE);
  fill_buff(buf[1],BUFSIZE);
  gmx_tx(dest,array(buf[1],BUFSIZE));
  gmx_rx(src,array(buf[0],BUFSIZE));
  gmx_tx_wait(dest);
  gmx_rx_wait(src);
  return check_buff(fp,"com_test",buf[0],BUFSIZE,1);
}

#ifdef USE_MAIN
#define	i860main_test	i860main
#endif

int i860main_test(int argc,char *argv[],FILE *stdlog,
                  int nnodes,int nodeid,int left,int right)
{
  int err;

  (void) fprintf(stdlog,"in i860main nnodes=%d, nodeid=%d, left=%d, right=%d, "
                 "bufsize=%d\n",nnodes,nodeid,left,right,BUFSIZE*sizeof(int));
  if (nnodes==1) test_connectivity(stdlog,sys_stat);
  err=com_test(stdlog,"left to right",left,right);
  err+=com_test(stdlog,"right to left",right,left);
  (void) fprintf(stdlog,"%d error(s)\n",err);
  fflush(stdlog);
  return 0;
}

/*
 * Debug routines end
 */

#ifdef DEBUG

static void print_stat(FILE *fp,char *title,t_stat *stat)
{
  (void) fprintf(fp,"%-35s (%s): cmd=%d, len=%6d, data=0x%x\n",
                 title,stat->name,peek_io_buf(IDT7130,stat->cmd),
                 stat->len,stat->data);
  fflush(fp);
}

#define ENTER(where,stat) print_stat(stdlog,"enter "#where,stat)
#define LEAVE(where,stat) print_stat(stdlog,"leave "#where,stat)
#else
#define	ENTER(where,stat)
#define	LEAVE(where,stat)
#endif

static void set_comm_led(int busy)
{
  if (busy)
    {
      if ((idt_busy++)==0) put_led(COMM_LED,1);
    }
  else
    {
      if ((--idt_busy)==0) put_led(COMM_LED,0);
    }
}

static int compare_tags(int tag1,int tag2)
{
  int i,count;

  count=0;
  for (i=-1; i!=0; i<<=1)
    {
      if ((tag1&1)==(tag2&1)) count++;
      tag1>>=1;
      tag2>>=1;      
    }
  return count;
}

void idtio_errstat(FILE *fp,char *msg,t_stat *stat,int dump)
{
  print_status(fp,idt_id,msg,sys_stat,stat);
  if (dump)
    {
      idt_dump(fp,"idtio_stat left",IDTLEFT,IRQBASE);
      idt_dump(fp,"idtio_stat right",IDTRIGHT,IRQBASE);
    }
  fflush(fp);
}

static int repair(FILE *fp,char *where,t_stat *stat,int tag)
{
  int out_tag,busy_tag,done_tag;

  busy_tag=compare_tags(tag,BUF_BUSY);
  done_tag=compare_tags(tag,BUF_DONE);
  if (busy_tag>done_tag)
    out_tag=BUF_BUSY;
  else
    if (done_tag>busy_tag)
      out_tag=BUF_DONE;
    else
      {
        out_tag=0;
        (void) fprintf(fp,"could not repair tag 0x%.2x in %s (%s)\n",tag,where,
                      stat->name);
        idtio_errstat(fp,where,stat,1);
        fatal_error(0,"received 0x%x as a tag on node %d",tag,idt_id);
      }
  (void) fprintf(fp,"repaired 0x%.2x in %s(%s) -> 0x%.2x\n",
                 tag,where,stat->name,out_tag);
  return out_tag;
}

static void buf_done(t_stat *stat)
{
  ENTER(buf_done,stat);
  set_comm_led(0);
  stat->data=NULL;
  stat->work=buf_idle;
  LEAVE(buf_done,stat);
}

static void buf_idle(t_stat *stat)
{
  ENTER(buf_idle,stat);
  LEAVE(buf_idle,stat);
}

static void buf_accept(t_stat *stat)
{
  int tag,repaired;
  static void buf_send(t_stat *stat);
  
  ENTER(buf_accept,stat);
  tag=peek_io_buf(IDT7130,stat->cmd);
  do
    {
      repaired=0;
      switch(tag)
        {
        case BUF_BUSY:
          stat->icnt++;
          break;
        case BUF_DONE:
          if (stat->len!=0) buf_send(stat); else buf_done(stat);
          break;
        default:
          tag=repair(stdlog,"buf_accept",stat,tag);
          repaired=1;
          break;
        }
    }
  while (repaired);
  LEAVE(buf_accept,stat);
}

static void buf_send(t_stat *stat)
{
  int len;
  
  ENTER(buf_send,stat);
  len=MIN(stat->bsize,stat->len);
  put_io_buf(IDT7130,stat->buf,stat->data,len);
  poke_io_buf(IDT7130,stat->cmd,BUF_BUSY);
  stat->bcnt++;
  stat->len-=len;
  stat->data+=len;
  stat->work=buf_accept;
  LEAVE(buf_send,stat);
}

static void buf_avail(t_stat *stat)
{
  int tag,len,repaired;
  
  ENTER(buf_avail,stat);
  tag=peek_io_buf(IDT7130,stat->cmd);
  do
    {
      repaired=0;
      switch (tag)
        {
        case BUF_DONE:
          stat->icnt++;
          break;
        case BUF_BUSY:
          len=MIN(stat->bsize,stat->len);
          get_io_buf(IDT7130,stat->buf,stat->data,len);
          poke_io_buf(IDT7130,stat->cmd,BUF_DONE);
          stat->bcnt++;
          stat->len-=len;
          stat->data+=len;
          if (stat->len==0) buf_done(stat);
          break;
        default:
          tag=repair(stdlog,"buf_avail",stat,tag);
          repaired=1;
          break;
        }
    }
  while (repaired);
  LEAVE(buf_avail,stat);
}

static void communicate(t_stat stat[IDTNR][TXRXNR])
{
  int i,j;
  
  for (i=0; i<IDTNR; i++)
    for (j=0; j<TXRXNR; j++)
      stat[i][j].work(&stat[i][j]);
}

static void chan_wait(t_stat stat[IDTNR][TXRXNR],int chan,int txrx)
{
  while (stat[chan][txrx].work!=buf_idle) communicate(stat);
}

static void fill_stat(t_stat *stat,void *bufptr,int bufsize,
                      void (*work)(t_stat *stat))
{
  ENTER(fill_stat,stat);
  stat->tcnt++;
  stat->len=bufsize;
  stat->data=bufptr;
  stat->work=work;
  set_comm_led(1);
  LEAVE(fill_stat,stat);
}

void idtio_tx(int chan,void *bufptr,int bufsize)
{
  ENTER(idtio_tx,&sys_stat[chan][TX]);
  chan_wait(sys_stat,chan,TX); /* Be sure that this channel is free for use */
  fill_stat(&sys_stat[chan][TX],bufptr,bufsize,buf_send);
  LEAVE(idtio_tx,&sys_stat[chan][TX]);
}

void idtio_rx(int chan,void *bufptr,int bufsize)
{
  ENTER(idtio_rx,&sys_stat[chan][RX]);
  chan_wait(sys_stat,chan,RX);
  fill_stat(&sys_stat[chan][RX],bufptr,bufsize,buf_avail);
  LEAVE(idtio_rx,&sys_stat[chan][RX]);
}

void idtio_tx_wait(int chan)
{
  ENTER(idtio_tx_wait,&sys_stat[chan][TX]);
  chan_wait(sys_stat,chan,TX);
  LEAVE(idtio_tx_wait,&sys_stat[chan][TX]);
}

void idtio_rx_wait(int chan)
{
  ENTER(idtio_rx_wait,&sys_stat[chan][RX]);
  chan_wait(sys_stat,chan,RX);
  LEAVE(idtio_rx_wait,&sys_stat[chan][RX]);
}

void idtio_txs(int chan,void *bufptr,int bufsize)
{
  ENTER(idtio_txs,&sys_stat[chan][RX]);
  idtio_tx(chan,bufptr,bufsize);
  idtio_tx_wait(chan);
  LEAVE(idtio_txs,&sys_stat[chan][RX]);
}

void idtio_rxs(int chan,void *bufptr,int bufsize)
{
  ENTER(idtio_rxs,&sys_stat[chan][RX]);
  idtio_rx(chan,bufptr,bufsize);
  idtio_rx_wait(chan);
  LEAVE(idtio_rxs,&sys_stat[chan][RX]);
}

static void init_stat(t_stat *stat,char *name,int bsize,int buf,int cmd)
{
  stat->name=name;
  stat->cmd=cmd;
  stat->buf=buf;
  stat->len=0;
  stat->bsize=bsize;
  stat->bcnt=0;
  stat->tcnt=0;
  stat->icnt=0;
  stat->data=NULL;
  stat->work=buf_idle;
  LEAVE(init_stat,stat);
}

static void ring_wait(int nodeid,int nnodes)
     /*
      * Wait for every node to execute this function, assume that nodeid
      * runs from 0 to nnodes - 1.
      */
{
  int i,newserver;
  char *sv,ch;

  if ((sv=getenv(SERVER_VERSION))==NULL) 
    newserver=0;
  else
    newserver=(strcmp(sv,USE_VERSION1)==0);
  if (newserver)
    syncall();
  else
    {
      ch=0;
      if (nodeid==0) 
        {
          for (i=1; i<nnodes; i++) rdlinda1in(i,-1,&ch);
          for (i=1; i<nnodes; i++) rdlinda1out(i+nnodes,1,&ch);
        }
      else
        {
          rdlinda1out(nodeid,1,&ch);
          rdlinda1in(nodeid+nnodes,-1,&ch);
        }
    }
}

void idtio_init(int nodeid,int nnodes)
{
  int i;

  idt_id=nodeid;
  init_stat(&sys_stat[LEFT] [TX],stat_names[LEFT][TX],
            BSIZE&~0x3,IDTLEFT,       IDTLEFT+2*BSIZE);
  init_stat(&sys_stat[LEFT] [RX],stat_names[LEFT][RX],
            BSIZE&~0x3,IDTLEFT+BSIZE, IDTLEFT+2*BSIZE+1);
  init_stat(&sys_stat[RIGHT][TX],stat_names[RIGHT][TX],
            BSIZE&~0x3,IDTRIGHT+BSIZE,IDTRIGHT+2*BSIZE+1);
  init_stat(&sys_stat[RIGHT][RX],stat_names[RIGHT][RX],
            BSIZE&~0x3,IDTRIGHT      ,IDTRIGHT+2*BSIZE);
  for (i=0; i<IRQBASE; i++) poke_io_buf(IDT7130,IDTLEFT+i,0xff);

#ifdef DEBUG
  print_status(stdlog,nodeid,"idtio_init (1)",sys_stat,NULL);
  ring_wait(nodeid,nnodes);
  print_status(stdlog,nodeid,"idtio_init (2)",sys_stat,NULL);
#endif
  ring_wait(nodeid,nnodes);

  for (i=0; i<IDTNR; i++) poke_io_buf(IDT7130,sys_stat[i][TX].cmd,BUF_DONE);
  ring_wait(nodeid,nnodes);
  for (i=0; i<IDTNR; i++) 
    {
      if (peek_io_buf(IDT7130,sys_stat[i][RX].cmd)!=BUF_DONE)
        {
          idtio_errstat(stdlog,"idtio_init",&sys_stat[i][RX],1);
          fatal_error(0,"idt %s not initialised on node %d",
                      sys_stat[i][RX].name,idt_id);
        }  
    }
#ifdef DEBUG
  print_status(stdlog,nodeid,"idtio_init(3)",sys_stat,NULL);
#endif
  ring_wait(nodeid,nnodes);
}

void idtio_stat(FILE *fp,char *msg)
{
  fprintf(fp,"idtio_stat message: %s\n",msg);
  fprintf(fp,"Idle Left  Send:    %d\n",sys_stat[LEFT][TX].icnt);
  fprintf(fp,"Idle Left  Receive: %d\n",sys_stat[LEFT][TX].icnt);
  fprintf(fp,"Idle Right Send:    %d\n",sys_stat[RIGHT][RX].icnt);
  fprintf(fp,"Idle Right Receive: %d\n",sys_stat[RIGHT][RX].icnt);
#ifdef DEBUG
  idtio_errstat(fp,msg,NULL,0);
#endif
}

void idt_left_right(int nnodes,int nodeid,int *left,int *right)
{
  *left=LEFT;
  *right=RIGHT;
}

void idt_reset_idle()
{
  sys_stat[LEFT][TX].icnt=0;
  sys_stat[LEFT][TX].icnt=0;
  sys_stat[RIGHT][RX].icnt=0;
  sys_stat[RIGHT][RX].icnt=0;
}

