#include <iostream>

class AbstractA {
public:
  AbstractA() { std::cout << "AbstractA default constructor" << std::endl; }
  AbstractA(std::string abstract_a_str, int abstract_a_int)
    : abstract_a_str_(abstract_a_str), abstract_a_int_(abstract_a_int)
  {
    std::cout << "AbstractA constructor" << std::endl;
  }
  // virtual destructor is needed to avoid memory leak
  // 基类的析构函数必须是虚函数！！！
  // 否则delete基类指针时，只会调用基类的析构函数，而不会调用派生类的析构函数，导致派生类的资源泄漏
  // 这里符合继承的关系，也就是父类无法调用子类的析构函数，所以需要将父类的析构函数设置为虚函数
  // 这样在delete基类指针时，会调用子类的析构函数
  virtual ~AbstractA() { std::cout << "AbstractA destructor" << std::endl; }
  // another way to define a pure virtual function -- default implementation
  void print() { std::cout << "Hello AbstractA!" << std::endl; }
  virtual void print_a() = 0;
  void print_abstract_a_int()
  {
    std::cout << "abstract_a_int_: " << abstract_a_int_ << std::endl;

    // using "this" to access the object's member variables
    std::cout << "abstract_a_int_: " << this->abstract_a_int_ << std::endl;
  }

  // return *this; 一定对应返回当前对象的引用
  // 因此可以继续对该对象进行操作
  AbstractA& set_abstract_a_int(int abstract_a_int_temp)
  {

    abstract_a_int_ = abstract_a_int_temp;

    // using "this" to access the object's member variables
    // this->abstract_a_int_ = abstract_a_int;

    return *this;
    // this是指向当前对象的指针
    // *this（解引用）是当前对象
    // 返回的是对象的引用
    // 这样可以实现链式调用
  }

  int get_a_int() { return abstract_a_int_; }

private:
  std::string abstract_a_str_;
  int abstract_a_int_;
};

class SubA : public AbstractA {
public:
  SubA() { std::cout << "SubA default constructor" << std::endl; }
  SubA(std::string abstract_a_str, int abstract_a_int, std::string sub_a_str,
    int sub_a_int)
    : AbstractA(abstract_a_str, abstract_a_int), sub_a_str_(sub_a_str),
      sub_a_int_(sub_a_int)
  {
    std::cout << "SubA constructor" << std::endl;
  }
  ~SubA() { std::cout << "SubA destructor" << std::endl; }
  void print() { std::cout << "Hello SubA!" << std::endl; }
  // the override can tell the compiler to check if the function is actually
  // overriding a virtual function they must have the same name and parameters
  // if it is not, the compiler will throw an error
  void print_a() override
  {
    std::cout << "Hello from SubA: " << sub_a_str_ << std::endl;
    // 基类的成员变量是private的，派生类无法直接访问
    // std::cout << "Hello from AbstractA: " << abstract_a_str_ << std::endl;
  }
  int get_a_int() { return sub_a_int_; }
  void print_different_a_int()
  {
    std::cout << AbstractA::get_a_int() << std::endl;
    std::cout << get_a_int() << std::endl;
  }

private:
  std::string sub_a_str_;
  int sub_a_int_;
};

int main()
{
  // initialize variables
  std::string abs_a_string = "abstract_a_string";
  int abs_a_integer = 1;
  std::string sub_a_string = "sub_a_string";
  int sub_a_integer = 2;

  // Error: cannot declare variable 'obj' to be of abstract type 'AbstractA'
  // AbstractA obj(abstract_a_str, abstract_a_int);

  // obj is an object of class SubA
  // stack allocation
  // 这里初始化是分配在栈上的，所以不需要delete
  SubA obj(abs_a_string, abs_a_integer, sub_a_string, sub_a_integer);

  obj.print();
  obj.print_a();
  obj.print_different_a_int();

  // 没有链式调用
  //
  // obj.set_abstract_a_int(20);
  // obj.print_abstract_a_int();
  // 链式调用（chain call）
  obj.set_abstract_a_int(20).print_abstract_a_int();

  std::cout << std::endl;

  // abstract_obj is a reference to obj
  // it's the AbstractA type,
  AbstractA& abstract_obj = obj;
  abstract_obj.print();
  abstract_obj.print_a();
  // error: 'class AbstractA' has no member named 'print_different_a_int'
  // abstract_obj.print_different_a_int();
  std::cout << std::endl;

  // sub_obj_ptr is a pointer to obj's address
  // it is A type
  SubA* sub_obj_ptr = &obj;
  sub_obj_ptr->print();
  sub_obj_ptr->print_a();
  std::cout << std::endl;

  // obj_abstract_ptr is a pointer to sub_obj_ptr
  // it is AbstractA type
  AbstractA* obj_abstract_ptr = sub_obj_ptr;
  obj_abstract_ptr->print();
  obj_abstract_ptr->print_a();
  std::cout << std::endl;

  // obj_abstract_ref is a reference to *sub_obj_ptr
  // it is AbstractA type
  AbstractA& obj_abstract_ref = *sub_obj_ptr;
  obj_abstract_ref.print();
  obj_abstract_ref.print_a();
  std::cout << std::endl;

  // heap allocation
  // 这里初始化是分配在堆上的，所以需要delete
  // 同时要注意基类的析构函数必须是虚函数，否则delete基类指针时，只会调用基类的析构函数，而不会调用派生类的析构函数
  AbstractA* heap_abstract_obj =
    new SubA(abs_a_string, abs_a_integer, sub_a_string, sub_a_integer);
  // AbstractA* heap_abstract_obj = new SubA();
  delete heap_abstract_obj;

  return 0;
}