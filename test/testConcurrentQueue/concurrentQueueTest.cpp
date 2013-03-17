#include "core.hpp"

using namespace lpt;

class Producer {
public:
   Producer(int id, concurrent_queue< lpt::ImageFrameGroup>& queue) : id(id), queue(queue), total(0)
   {
      image_type = cv::Mat(cv::Size(640,480), CV_8UC1);
   }
   virtual void operator()() {
      cout << "Producer " << id << " starting" << endl;

      boost::posix_time::microseconds worktime(1);
      for(int i = 0; i < 500; i++) {
         lpt::ImageFrameGroup group(6);
         for (int j = 0; j < 6; j++) {
            group[j].image = image_type.clone();
            group[j].frame_index = i;
         }
         if (queue.push(group))
         {
            ++total;
         }
         boost::this_thread::sleep(worktime);
      }
      cout << "Producer " << id << " done" << endl;
   }

   unsigned int numProduced(){ return total;}

private:
   int id;
   unsigned int total;
   concurrent_queue<lpt::ImageFrameGroup>& queue;
   cv::Mat image_type;
};

class Consumer {
public: 
   Consumer(int id, concurrent_queue< lpt::ImageFrameGroup>& queue) : id(id), queue(queue), total(0) {}
   virtual void operator()() 
   {
      cout << "Consumer " << id << " starting " << endl;
      boost::posix_time::milliseconds worktime(2);
      try
      {
         while ( 1 ) 
         {
            lpt::ImageFrameGroup group;	
            queue.wait_and_pop(group);
            ++total;
            boost::this_thread::sleep(worktime);
            boost::this_thread::interruption_point();
         }
      }
      catch (boost::thread_interrupted)
      {
         cout << "Consumer " << id << " interupted " << endl;
      }
   }

   unsigned int numConsumed(){ return total;}
private:
   int id;
   unsigned int total;
   concurrent_queue<lpt::ImageFrameGroup>& queue;
};

int main() {

   concurrent_queue< lpt::ImageFrameGroup >  buf;
   buf.setCapacity(500);

   int num_producers = 5;
   int num_consumers = 3;

   boost::thread_group producer_group;
   boost::thread_group consumer_group;

   vector<Producer> producers;
   vector<Consumer> consumers;

   for (int i = 0; i < num_producers; ++i) {
      producers.push_back(Producer (i, buf));
      producer_group.create_thread(producers.back());
   }

   for (int i = 0; i < num_consumers; ++i) {
      consumers.push_back(Consumer(i,buf));
      consumer_group.create_thread(consumers.back());
   }

   producer_group.join_all();

   boost::posix_time::milliseconds worktime(1000);
   boost::this_thread::sleep(worktime);

   cout << "Main thread interupting consumers" << endl;
   consumer_group.interrupt_all();
   consumer_group.join_all();
   unsigned int num_produced = 0;
   for (int i = 0; i < producers.size(); ++i) {
      num_produced += producers[i].numProduced();
   }

   unsigned int num_consumed = 0;
   for (int i = 0; i < num_consumers; ++i) {
      num_consumed += consumers[i].numConsumed();
   }
   bool passed = num_produced == num_consumed + buf.size() ? true : false;
   cout << "Test " << (passed ? "PASSED" : "FAILED") << endl;
   return 0;
}


