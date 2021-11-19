#include <atomic>
#include <iostream>
#include <thread>

/* bool atomic<T>::compare_exchange_weak(T& expected, T desired)
 * 現在の値が expected に等しければ desired を書いて true を返す
 * そうでなければ expected を現在の値に書き換えて false を返す
 * これはexchangeが失敗する可能性があるため、loopが必要となる。
 */

/* bool atomic<T>::compare_exchange_strong(T& expected, T desired)
 * 現在の値が expected に等しければ desired を書いて true を返す
 * そうでなければ expected を現在の値に書き換えて false を返す
 * exchangeは必ず成功するためloopを使わない場合に使ってね。
 */

int main()
{
  std::atomic<int> val(0);
  
  std::thread th1([&]() {
                    for(int i=0; i<10000000; i++){
                      int expected = val.load();
                      int desired;
                      do {
                        desired = expected + 1;
                      } while ( !val.compare_exchange_weak(expected, desired));
                    }
                  });
  std::thread th2([&]() {
                    for (int i = 0; i < 10000000; i++) {
                      int expected = val.load();
                      int desired;
                      do {
                        desired = expected - 1;
                      } while (!val.compare_exchange_weak(expected, desired));
                    }
                  });
  std::thread th3([&]() {
                    for (int i = 0; i < 10000000; i++) {
                      int expected = val.load();
                      int desired;
                      do {
                        desired = expected - 1;
                      } while (!val.compare_exchange_weak(expected, desired));
                    }
                  });
  th1.join();
  th2.join();
  th3.join();
  std::cout << val << std::endl;
}
