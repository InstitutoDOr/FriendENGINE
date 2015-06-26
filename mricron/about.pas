unit about;

interface

uses
{$IFDEF FPC}LResources,{$ELSE} ShellAPI, {$ENDIF}
{$IFNDEF Unix} Windows,{$ENDIF}
   SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, ExtCtrls, define_types;

type

  { TAboutForm }

  TAboutForm = class(TForm)
    HomepageLabel: TLabel;
    Panel2: TPanel;
    Panel1: TPanel;
    ThreadLabel: TLabel;
    procedure FormCreate(Sender: TObject);
    procedure HomePageClick(Sender: TObject);
    procedure Panel1Click(Sender: TObject);
    procedure Panel1DragDrop(Sender, Source: TObject; X, Y: Integer);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  AboutForm: TAboutForm;

implementation
{$IFNDEF FPC}
{$R *.DFM}
{$ENDIF}

procedure TAboutForm.FormCreate(Sender: TObject);
begin
      HomepageLabel.caption := 'www.mricro.com :: '+kMRIcronVers ;
end;

procedure TAboutForm.HomePageClick(Sender: TObject);
begin
{$IFDEF FPC}
{$ELSE}
   ShellExecute (0, Nil, 'http://www.mricro.com', Nil, Nil, SW_ShowDefault);
{$ENDIF}
end;

procedure TAboutForm.Panel1Click(Sender: TObject);
begin

end;

procedure TAboutForm.Panel1DragDrop(Sender, Source: TObject; X, Y: Integer);
begin
//showmessage('x');
end;

{$IFDEF FPC}
initialization
  {$I about.lrs}
{$ENDIF}

end.
