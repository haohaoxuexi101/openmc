#include <iostream>

class AbstractA {
public:
  AbstractA(std::string a_str, int a_int) : a_str_(a_str), a_int_(a_int) {}
  void print() { std::cout << "Hello from AbstractA " << a_str_ << std::endl; }
  virtual void print_a() = 0;
  void print_a_int() { std::cout << "a_int_: " << a_int_ << std::endl; }
  // using "this" to access the object's member variables
  void print_a_int_this()
  {
    std::cout << "a_int_: " << this->a_int_ << std::endl;
  }
  AbstractA& set_a_int(int a_int)
  {
    // this->a_int_ = a_int;
    a_int_ = a_int;
    return *this;
  }
  void set_a_int_this(int a_int)
  {
    // this->a_int_ = a_int;
    a_int_ = a_int;
  }
  int get_a_int() { return a_int_; }

private:
  std::string a_str_;
  int a_int_;
};

class A : public AbstractA {
public:
  A(std::string a_str, int a_int, std::string new_a_str, int new_a_int)
    : AbstractA(a_str, a_int), new_a_str_(new_a_str), new_a_int_(new_a_int),
      a_int_(new_a_int)
  {}
  void print() { std::cout << "Hello from A" << std::endl; }
  // the override can tell the compiler to check if the function is actually
  // overriding a virtual function they must have the same name and parameters
  // if it is not, the compiler will throw an error
  void print_a() override { std::cout << "Hello from A" << std::endl; }
  int get_a_int() { return a_int_; }
  void print_a_dif()
  {
    std::cout << AbstractA::get_a_int() << std::endl;
    std::cout << get_a_int() << std::endl;
  }

private:
  std::string new_a_str_;
  int new_a_int_;
  int a_int_;
};

int main()
{
  // initialize variables
  std::string a_str = "a_string";
  int a_int = 1;

  // Error: cannot declare variable 'obj' to be of abstract type 'AbstractA'
  // AbstractA obj(a_str, a_int);

  // obj is an object of class A
  A obj(a_str, a_int, a_str, 2);

  obj.print();
  obj.print_a();
  obj.print_a_dif();

  // obj.set_a_int(2);
  // obj.set_a_int_this(3);
  // obj.print_a_int();
  // obj.print_a_int_this();
  std::cout << std::endl;

  // obj_abstract is a reference to obj
  // it's the AbstractA type,
  AbstractA& obj_abstract = obj;
  obj_abstract.print();
  obj_abstract.print_a();
  // error: 'class AbstractA' has no member named 'print_a_diff'
  // obj_abstract.print_a_dif();
  std::cout << std::endl;

  // obj_ptr is a pointer to obj's address
  // it is A type
  A* obj_ptr = &obj;
  obj_ptr->print();
  obj_ptr->print_a();
  std::cout << std::endl;

  // obj_abstract_ptr is a pointer to obj_ptr
  // it is AbstractA type
  AbstractA* obj_abstract_ptr = obj_ptr;
  obj_abstract_ptr->print();
  obj_abstract_ptr->print_a();
  std::cout << std::endl;

  // obj_abstract_ref is a reference to *obj_ptr
  // it is AbstractA type
  AbstractA& obj_abstract_ref = *obj_ptr;
  obj_abstract_ref.print();
  obj_abstract_ref.print_a();
  std::cout << std::endl;

  return 0;
}